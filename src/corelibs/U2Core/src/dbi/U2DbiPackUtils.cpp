/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2017 UniPro <ugene@unipro.ru>
 * http://ugene.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

#include <U2Core/DatatypeSerializeUtils.h>
#include <U2Core/U2Feature.h>
#include <U2Core/U2Region.h>
#include <U2Core/U2OpStatusUtils.h>
#include <U2Core/U2SafePoints.h>

#include "U2DbiPackUtils.h"

namespace U2 {

const QByteArray U2DbiPackUtils::VERSION("0");
const char U2DbiPackUtils::SEP = '\t';
const char U2DbiPackUtils::SECOND_SEP = 11;
const char U2DbiPackUtils::THIRD_SEP = 12;

QByteArray U2DbiPackUtils::packGaps(const QList<U2MsaGap> &gaps) {
    QByteArray result;
    foreach (const U2MsaGap &gap, gaps) {
        if (!result.isEmpty()) {
            result += ";";
        }
        result += QByteArray::number(gap.offset);
        result += ",";
        result += QByteArray::number(gap.gap);
    }
    return "\"" + result + "\"";
}

bool U2DbiPackUtils::unpackGaps(const QByteArray &str, QList<U2MsaGap> &gaps) {
    CHECK(str.startsWith('\"') && str.endsWith('\"'), false);
    QByteArray gapsStr = str.mid(1, str.length() - 2);
    if (gapsStr.isEmpty()) {
        return true;
    }

    QList<QByteArray> tokens = gapsStr.split(';');
    foreach (const QByteArray &t, tokens) {
        QList<QByteArray> gapTokens = t.split(',');
        CHECK(2 == gapTokens.size(), false);
        bool ok = false;
        U2MsaGap gap;
        gap.offset = gapTokens[0].toLongLong(&ok);
        CHECK(ok, false);
        gap.gap = gapTokens[1].toLongLong(&ok);
        CHECK(ok, false);
        gaps << gap;
    }
    return true;
}

QByteArray U2DbiPackUtils::packGapDetails(qint64 rowId, const QList<U2MsaGap> &oldGaps, const QList<U2MsaGap> &newGaps) {
    QByteArray result = VERSION;
    result += SEP;
    result += QByteArray::number(rowId);
    result += SEP;
    result += packGaps(oldGaps);
    result += SEP;
    result += packGaps(newGaps);
    return result;
}

bool U2DbiPackUtils::unpackGapDetails(const QByteArray &modDetails, qint64 &rowId, QList<U2MsaGap> &oldGaps, QList<U2MsaGap> &newGaps) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(4 == tokens.size(), QString("Invalid gap modDetails string '%1'").arg(QString(modDetails)), false);
    { // version
        SAFE_POINT(VERSION == tokens[0], QString("Invalid modDetails version '%1'").arg(tokens[0].data()), false);
    }
    { // rowId
        bool ok = false;
        rowId = tokens[1].toLongLong(&ok);
        SAFE_POINT(ok, QString("Invalid gap modDetails rowId '%1'").arg(tokens[1].data()), false);
    }
    { // oldGaps
        bool ok = unpackGaps(tokens[2], oldGaps);
        SAFE_POINT(ok, QString("Invalid gap string '%1'").arg(tokens[2].data()), false);
    }
    { // newGaps
        bool ok = unpackGaps(tokens[3], newGaps);
        SAFE_POINT(ok, QString("Invalid gap string '%1'").arg(tokens[3].data()), false);
    }
    return true;
}

QByteArray U2DbiPackUtils::packGapDetails(qint64 rowId, const U2DataId &relatedObjectId, const QList<U2MsaGap> &oldGaps, const QList<U2MsaGap> &newGaps) {
    QByteArray result = VERSION;
    result += SEP;
    result += QByteArray::number(rowId);
    result += SEP;
    result += relatedObjectId.toHex();
    result += SEP;
    result += packGaps(oldGaps);
    result += SEP;
    result += packGaps(newGaps);
    return result;
}

bool U2DbiPackUtils::unpackGapDetails(const QByteArray &modDetails, qint64 &rowId, U2DataId &relatedObjectId, QList<U2MsaGap> &oldGaps, QList<U2MsaGap> &newGaps) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(5 == tokens.size(), QString("Invalid gap modDetails string '%1'").arg(QString(modDetails)), false);
    { // version
        SAFE_POINT(VERSION == tokens[0], QString("Invalid modDetails version '%1'").arg(tokens[0].data()), false);
    }
    { // rowId
        bool ok = false;
        rowId = tokens[1].toLongLong(&ok);
        SAFE_POINT(ok, QString("Invalid gap modDetails rowId '%1'").arg(tokens[1].data()), false);
    }
    { // relatedObjectId
        relatedObjectId = QByteArray::fromHex(tokens[2]);
    }
    { // oldGaps
        bool ok = unpackGaps(tokens[3], oldGaps);
        SAFE_POINT(ok, QString("Invalid gap string '%1'").arg(tokens[3].data()), false);
    }
    { // newGaps
        bool ok = unpackGaps(tokens[4], newGaps);
        SAFE_POINT(ok, QString("Invalid gap string '%1'").arg(tokens[4].data()), false);
    }
    return true;
}

QByteArray U2DbiPackUtils::packRowOrder(const QList<qint64>& rowIds) {
    QByteArray result;
    foreach (qint64 rowId, rowIds) {
        if (!result.isEmpty()) {
            result += ",";
        }
        result += QByteArray::number(rowId);
    }
    return "\"" + result + "\"";
}

bool U2DbiPackUtils::unpackRowOrder(const QByteArray& str, QList<qint64>& rowsIds) {
    CHECK(str.startsWith('\"') && str.endsWith('\"'), false);
    QByteArray orderStr = str.mid(1, str.length() - 2);
    if (orderStr.isEmpty()) {
        return true;
    }

    QList<QByteArray> tokens = orderStr.split(',');
    foreach (const QByteArray &t, tokens) {
        bool ok = false;
        rowsIds << t.toLongLong(&ok);
        CHECK(ok, false);
    }
    return true;
}

QByteArray U2DbiPackUtils::packRowOrderDetails(const QList<qint64>& oldOrder, const QList<qint64>& newOrder) {
    QByteArray result = VERSION;
    result += SEP;
    result += packRowOrder(oldOrder);
    result += SEP;
    result += packRowOrder(newOrder);
    return result;
}

bool U2DbiPackUtils::unpackRowOrderDetails(const QByteArray &modDetails, QList<qint64>& oldOrder, QList<qint64>& newOrder) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(3 == tokens.size(), QString("Invalid rows order modDetails string '%1'").arg(QString(modDetails)), false);
    { // version
        SAFE_POINT(VERSION == tokens[0], QString("Invalid modDetails version '%1'").arg(tokens[0].data()), false);
    }
    { // oldOrder
        bool ok = unpackRowOrder(tokens[1], oldOrder);
        SAFE_POINT(ok, QString("Invalid rows order string '%1'").arg(tokens[1].data()), false);
    }
    { // newGaps
        bool ok = unpackRowOrder(tokens[2], newOrder);
        SAFE_POINT(ok, QString("Invalid rows order string '%1'").arg(tokens[2].data()), false);
    }

    return true;
}

QByteArray U2DbiPackUtils::packRowNameDetails(qint64 rowId, const QString &oldName, const QString &newName) {
    QByteArray result = VERSION;
    result += SEP;
    result += QByteArray::number(rowId);
    result += SEP;
    result += oldName.toLatin1();
    result += SEP;
    result += newName.toLatin1();
    return result;
}

bool U2DbiPackUtils::unpackRowNameDetails(const QByteArray &modDetails, qint64 &rowId, QString &oldName, QString &newName) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(4 == tokens.size(), QString("Invalid row name modDetails string '%1'").arg(QString(modDetails)), false);
    { // version
        SAFE_POINT(VERSION == tokens[0], QString("Invalid modDetails version '%1'").arg(tokens[0].data()), false);
    }
    { // rowId
        bool ok = false;
        rowId = tokens[1].toLongLong(&ok);
        SAFE_POINT(ok, QString("Invalid row name modDetails rowId '%1'").arg(tokens[1].data()), false);
    }
    { // oldName
        oldName = tokens[2];
    }
    { // newName
        newName = tokens[3];
    }
    return true;
}

QByteArray U2DbiPackUtils::packRow(qint64 posInMsa, const U2MsaRow& row) {
    QByteArray result = VERSION;
    result += SEP;
    result += QByteArray::number(posInMsa);
    result += SEP;
    result += QByteArray::number(row.rowId);
    result += SEP;
    result += row.sequenceId.toHex();
    result += SEP;
    result += QByteArray::number(row.gstart);
    result += SEP;
    result += QByteArray::number(row.gend);
    result += SEP;
    result += packGaps(row.gaps);
    return result;
}

bool U2DbiPackUtils::unpackRow(const QByteArray &modDetails, qint64& posInMsa, U2MsaRow& row) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(7 == tokens.size(), QString("Invalid added row modDetails string '%1'").arg(QString(modDetails)), false);
    { // version
        SAFE_POINT(VERSION == tokens[0], QString("Invalid modDetails version '%1'").arg(tokens[0].data()), false);
    }
    { // posInMsa
        bool ok = false;
        posInMsa = tokens[1].toLongLong(&ok);
        SAFE_POINT(ok, QString("Invalid added row modDetails posInMsa '%1'").arg(tokens[1].data()), false);
    }
    { // rowId
        bool ok = false;
        row.rowId = tokens[2].toLongLong(&ok);
        SAFE_POINT(ok, QString("Invalid added row modDetails rowId '%1'").arg(tokens[2].data()), false);
    }
    { // sequenceId
        row.sequenceId = QByteArray::fromHex(tokens[3]);
    }
    { // gstart
        bool ok = false;
        row.gstart = tokens[4].toLongLong(&ok);
        SAFE_POINT(ok, QString("Invalid added row modDetails gstart '%1'").arg(tokens[4].data()), false);
    }
    { // gend
        bool ok = false;
        row.gend = tokens[5].toLongLong(&ok);
        SAFE_POINT(ok, QString("Invalid added row modDetails gend '%1'").arg(tokens[5].data()), false);
    }
    { // gaps
        bool ok = unpackGaps(tokens[6], row.gaps);
        SAFE_POINT(ok, QString("Invalid added row modDetails gaps '%1'").arg(tokens[6].data()), false);
    }
    return true;
}

QByteArray U2DbiPackUtils::packRow(qint64 posInMca, const U2McaRow &row) {
    QByteArray result = VERSION;
    result += SEP;
    result += QByteArray::number(posInMca);
    result += SEP;
    result += QByteArray::number(row.rowId);
    result += SEP;
    result += row.chromatogramId.toHex();
    result += SEP;
    result += row.sequenceId.toHex();
    result += SEP;
    result += QByteArray::number(row.gstart);
    result += SEP;
    result += QByteArray::number(row.gend);
    result += SEP;
    result += packGaps(row.gaps);
    return result;
}

bool U2DbiPackUtils::unpackRow(const QByteArray &modDetails, qint64 &posInMca, U2McaRow &row) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(8 == tokens.size(), QString("Invalid added row modDetails string '%1'").arg(QString(modDetails)), false);
    { // version
        SAFE_POINT(VERSION == tokens[0], QString("Invalid modDetails version '%1'").arg(tokens[0].data()), false);
    }
    { // posInMsa
        bool ok = false;
        posInMca = tokens[1].toLongLong(&ok);
        SAFE_POINT(ok, QString("Invalid added row modDetails posInMsa '%1'").arg(tokens[1].data()), false);
    }
    { // rowId
        bool ok = false;
        row.rowId = tokens[2].toLongLong(&ok);
        SAFE_POINT(ok, QString("Invalid added row modDetails rowId '%1'").arg(tokens[2].data()), false);
    }
    { // chromatogramId
        row.chromatogramId = QByteArray::fromHex(tokens[3]);
    }
    { // sequenceId
        row.sequenceId = QByteArray::fromHex(tokens[4]);
    }
    { // gstart
        bool ok = false;
        row.gstart = tokens[5].toLongLong(&ok);
        SAFE_POINT(ok, QString("Invalid added row modDetails gstart '%1'").arg(tokens[5].data()), false);
    }
    { // gend
        bool ok = false;
        row.gend = tokens[6].toLongLong(&ok);
        SAFE_POINT(ok, QString("Invalid added row modDetails gend '%1'").arg(tokens[6].data()), false);
    }
    { // sequence gaps
        bool ok = unpackGaps(tokens[7], row.gaps);
        SAFE_POINT(ok, QString("Invalid added row modDetails gaps '%1'").arg(tokens[7].data()), false);
    }
    return true;
}

QByteArray U2DbiPackUtils::packRowInfo(const U2MsaRow &row) {
    QByteArray result;
    result += QByteArray::number(row.rowId);
    result += SECOND_SEP;
    result += row.sequenceId.toHex();
    result += SECOND_SEP;
    result += QByteArray::number(row.gstart);
    result += SECOND_SEP;
    result += QByteArray::number(row.gend);
    result += SECOND_SEP;
    result += QByteArray::number(row.length);
    return result;
}

bool U2DbiPackUtils::unpackRowInfo(const QByteArray &str, U2MsaRow& row) {
    QList<QByteArray> tokens = str.split(SECOND_SEP);
    CHECK(5 == tokens.count(), false);

    bool ok = false;

    row.rowId = tokens[0].toLongLong(&ok);
    CHECK(ok, false);
    row.sequenceId = QByteArray::fromHex(tokens[1]);
    row.gstart = tokens[2].toLongLong(&ok);
    CHECK(ok, false);
    row.gend = tokens[3].toLongLong(&ok);
    CHECK(ok, false);
    row.length = tokens[4].toLongLong(&ok);
    CHECK(ok, false);

    return true;
}

QByteArray U2DbiPackUtils::packRowInfo(const U2McaRow &row) {
    QByteArray result;
    result += QByteArray::number(row.rowId);
    result += SECOND_SEP;
    result += row.chromatogramId.toHex();
    result += SECOND_SEP;
    result += row.sequenceId.toHex();
    result += SECOND_SEP;
    result += QByteArray::number(row.gstart);
    result += SECOND_SEP;
    result += QByteArray::number(row.gend);
    result += SECOND_SEP;
    result += QByteArray::number(row.length);
    return result;
}

bool U2DbiPackUtils::unpackRowInfo(const QByteArray &str, U2McaRow &row) {
    QList<QByteArray> tokens = str.split(SECOND_SEP);
    CHECK(5 == tokens.count(), false);

    bool ok = false;

    row.rowId = tokens[0].toLongLong(&ok);
    CHECK(ok, false);
    row.chromatogramId = QByteArray::fromHex(tokens[1]);
    row.sequenceId = QByteArray::fromHex(tokens[2]);
    row.gstart = tokens[3].toLongLong(&ok);
    CHECK(ok, false);
    row.gend = tokens[4].toLongLong(&ok);
    CHECK(ok, false);
    row.length = tokens[5].toLongLong(&ok);
    CHECK(ok, false);

    return true;
}

QByteArray U2DbiPackUtils::packRowInfoDetails(const U2MsaRow &oldRow, const U2MsaRow &newRow) {
    QByteArray result = VERSION;
    result += SEP;
    result += packRowInfo(oldRow);
    result += SEP;
    result += packRowInfo(newRow);
    return result;
}

bool U2DbiPackUtils::unpackRowInfoDetails(const QByteArray &modDetails, U2MsaRow &oldRow, U2MsaRow &newRow) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(3 == tokens.count(), QString("Invalid modDetails '%1'!").arg(QString(modDetails)), false);
    SAFE_POINT(VERSION == tokens[0], QString("Invalid modDetails version '%1'").arg(QString(tokens[0])), false);

    bool ok = false;
    ok = unpackRowInfo(tokens[1], oldRow);
    CHECK(ok, false);
    ok = unpackRowInfo(tokens[2], newRow);
    CHECK(ok, false);

    return true;
}

QByteArray U2DbiPackUtils::packRowInfoDetails(const U2McaRow &oldRow, const U2McaRow &newRow) {
    QByteArray result = VERSION;
    result += SEP;
    result += packRowInfo(oldRow);
    result += SEP;
    result += packRowInfo(newRow);
    return result;
}

bool U2DbiPackUtils::unpackRowInfoDetails(const QByteArray &modDetails, U2McaRow &oldRow, U2McaRow &newRow) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(3 == tokens.count(), QString("Invalid modDetails '%1'").arg(QString(modDetails)), false);
    SAFE_POINT(VERSION == tokens[0], QString("Invalid modDetails version '%1'").arg(QString(tokens[0])), false);

    bool ok = false;
    ok = unpackRowInfo(tokens[1], oldRow);
    CHECK(ok, false);
    ok = unpackRowInfo(tokens[2], newRow);
    CHECK(ok, false);

    return true;
}

QByteArray U2DbiPackUtils::packAlphabetDetails(const U2AlphabetId &oldAlphabet, const U2AlphabetId &newAlphabet) {
    QByteArray result = VERSION;
    result += SEP;
    result += oldAlphabet.id.toLatin1();
    result += SEP;
    result += newAlphabet.id.toLatin1();
    return result;
}

bool U2DbiPackUtils::unpackAlphabetDetails(const QByteArray &modDetails, U2AlphabetId &oldAlphabet, U2AlphabetId &newAlphabet) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(3 == tokens.count(), QString("Invalid modDetails '%1'!").arg(QString(modDetails)), false);
    SAFE_POINT(VERSION == tokens[0], QString("Invalid modDetails version '%1'").arg(QString(tokens[0])), false);

    oldAlphabet = QString(tokens[1]);
    newAlphabet = QString(tokens[2]);

    if (!oldAlphabet.isValid() || !newAlphabet.isValid()) {
        return false;
    }

    return true;
}

QByteArray U2DbiPackUtils::packObjectNameDetails(const QString &oldName, const QString &newName) {
    QByteArray result = VERSION;
    result += SEP;
    result += oldName.toUtf8();
    result += SEP;
    result += newName.toUtf8();
    return result;
}

bool U2DbiPackUtils::unpackObjectNameDetails(const QByteArray &modDetails, QString &oldName, QString &newName) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(3 == tokens.count(), "Invalid modDetails!", false);
    SAFE_POINT(VERSION == tokens[0], "Invalid modDetails version!", false);
    SAFE_POINT(!QString(tokens[1]).isEmpty(), "Invalid modDetails!", false);
    SAFE_POINT(!QString(tokens[2]).isEmpty(), "Invalid modDetails!", false);

    oldName = QString::fromUtf8(tokens[1]);
    newName = QString::fromUtf8(tokens[2]);
    return true;
}

QByteArray U2DbiPackUtils::packRows(const QList<qint64> &posInMsa, const QList<U2MsaRow> &rows) {
    SAFE_POINT(posInMsa.size() == rows.size(), "Different lists sizes", "");
    QByteArray result = VERSION;
    QList<qint64>::ConstIterator pi = posInMsa.begin();
    QList<U2MsaRow>::ConstIterator ri = rows.begin();
    for (; ri != rows.end(); ri++, pi++) {
        result += SECOND_SEP + packRow(*pi, *ri);
    }
    return result;
}

bool U2DbiPackUtils::unpackRows(const QByteArray &modDetails, QList<qint64> &posInMsa, QList<U2MsaRow> &rows) {
    QList<QByteArray> tokens = modDetails.split(SECOND_SEP);
    SAFE_POINT(tokens.count() > 0, QString("Invalid modDetails '%1'!").arg(QString(modDetails)), false);
    QByteArray modDetailsVersion = tokens.takeFirst();
    SAFE_POINT(VERSION == modDetailsVersion, QString("Invalid modDetails version '%1'").arg(QString(modDetailsVersion)), false);
    foreach (const QByteArray &token, tokens) {
        qint64 pos = 0;
        U2MsaRow row;
        bool ok = unpackRow(token, pos, row);
        CHECK(ok, false);
        posInMsa << pos;
        rows << row;
    }
    return true;
}

QByteArray U2DbiPackUtils::packRows(const QList<qint64> &posInMca, const QList<U2McaRow> &rows) {
    SAFE_POINT(posInMca.size() == rows.size(), "Different lists sizes", "");
    QByteArray result = VERSION;
    QList<qint64>::ConstIterator pi = posInMca.begin();
    QList<U2McaRow>::ConstIterator ri = rows.begin();
    for (; ri != rows.end(); ri++, pi++) {
        result += SECOND_SEP + packRow(*pi, *ri);
    }
    return result;
}

bool U2DbiPackUtils::unpackRows(const QByteArray &modDetails, QList<qint64> &posInMca, QList<U2McaRow> &rows) {
    QList<QByteArray> tokens = modDetails.split(SECOND_SEP);
    SAFE_POINT(tokens.count() > 0, QString("Invalid modDetails '%1'").arg(QString(modDetails)), false);
    QByteArray modDetailsVersion = tokens.takeFirst();
    SAFE_POINT(VERSION == modDetailsVersion, QString("Invalid modDetails version '%1'").arg(QString(modDetailsVersion)), false);
    foreach (const QByteArray &token, tokens) {
        qint64 pos = 0;
        U2McaRow row;
        bool ok = unpackRow(token, pos, row);
        CHECK(ok, false);
        posInMca << pos;
        rows << row;
    }
    return true;
}

QByteArray U2DbiPackUtils::packSequenceDataDetails(const U2Region &replacedRegion, const QByteArray &oldData,
                                              const QByteArray &newData, const QVariantMap &hints) {
    // replacedRegion length is used only for check, it is not stored
    SAFE_POINT(replacedRegion.length >= oldData.length(), "oldData length does not match to the region length.", QByteArray());
    QByteArray result = VERSION;
    result += SEP;
    result += QByteArray::number(replacedRegion.startPos);
    result += SEP;
    result += oldData;
    result += SEP;
    result += newData;
    result += SEP;
    result += packSequenceDataHints(hints);
    return result;
}

bool U2DbiPackUtils::unpackSequenceDataDetails(const QByteArray &modDetails, U2Region &replacedRegion, QByteArray &oldData,
                                          QByteArray &newData, QVariantMap &hints) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(tokens.count() == 5, QString("Invalid modDetails '%1'!").arg(QString(modDetails)), false);
    SAFE_POINT(VERSION == tokens[0], QString("Invalid modDetails version '%1'").arg(QString(tokens[0])), false);

    SAFE_POINT(!QString(tokens[1]).isEmpty(), "Invalid modDetails!", false);

    bool ok = false;
    replacedRegion = U2Region(tokens[1].toLongLong(&ok), tokens[2].length());
    CHECK(ok, false);

    oldData = tokens[2];
    newData = tokens[3];
    ok = unpackSequenceDataHints(tokens[4], hints);
    CHECK(ok, false);
    return true;
}

QByteArray U2DbiPackUtils::packChromatogramData(const DNAChromatogram &chromatogram) {
    return DNAChromatogramSerializer::serialize(chromatogram).toHex();
}

bool U2DbiPackUtils::unpackChromatogramData(const QByteArray &modDetails, DNAChromatogram &chromatogram) {
    U2OpStatusImpl os;
    chromatogram = DNAChromatogramSerializer::deserialize(QByteArray::fromHex(modDetails), os);
    return !os.hasError();
}

QByteArray U2DbiPackUtils::packChromatogramDetails(const DNAChromatogram &oldChromatogram, const DNAChromatogram &newChromatogram) {
    QByteArray result = VERSION;
    result += SEP;
    result += packChromatogramData(oldChromatogram);
    result += SEP;
    result += packChromatogramData(newChromatogram);
    return result;
}

bool U2DbiPackUtils::unpackChromatogramDetails(const QByteArray &modDetails, DNAChromatogram &oldChromatogram, DNAChromatogram &newChromatogram) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(tokens.count() == 3, QString("Invalid modDetails '%1'").arg(QString(modDetails)), false);
    SAFE_POINT(VERSION == tokens[0], QString("Invalid modDetails version '%1'").arg(QString(tokens[0])), false);

    bool ok = false;
    ok = unpackChromatogramData(tokens[1], oldChromatogram);
    CHECK(ok, false);

    ok = unpackChromatogramData(tokens[2], newChromatogram);
    CHECK(ok, false);

    return true;
}

QByteArray U2DbiPackUtils::packSequenceDataHints(const QVariantMap &hints) {
    QByteArray result;
    foreach (QString key, hints.keys()) {
        if (false == result.isEmpty()) {
            result += ";";
        }
        result += key + "," + hints[key].toByteArray();
    }
    return "\"" + result + "\"";
}

bool U2DbiPackUtils::unpackSequenceDataHints(const QByteArray &str, QVariantMap &hints) {
    CHECK(str.startsWith('\"') && str.endsWith('\"'), false);
    QByteArray hintsStr = str.mid(1, str.length() - 2);
    if (hintsStr.isEmpty()) {
        return true;
    }

    QList<QByteArray> tokens = hintsStr.split(';');
    foreach (const QByteArray &t, tokens) {
        QList<QByteArray> hintTokens = t.split(',');
        CHECK(2 == hintTokens.size(), false);
        hints.insert(QString(hintTokens[0]), QVariant(hintTokens[1]));
    }
    return true;
}

QByteArray U2DbiPackUtils::packAlignmentLength(const qint64 oldLen, const qint64 newLen) {
    QByteArray result;
    result += QByteArray::number(oldLen);
    result += SEP;
    result += QByteArray::number(newLen);
    return result;
}

bool U2DbiPackUtils::unpackAlignmentLength(const QByteArray &modDetails, qint64 &oldLen, qint64 &newLen) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(tokens.count() == 2, QString("Invalid modDetails '%1'!").arg(QString(modDetails)), false);

    bool ok = false;
    oldLen = tokens.first().toInt(&ok);
    CHECK(ok, false);
    newLen = tokens.last().toInt(&ok);
    CHECK(ok, false);

    return true;
}

QByteArray U2DbiPackUtils::packUdr(const QByteArray& oldData, const QByteArray& newData) {
    QByteArray result;
    result += oldData.toHex();
    result += SEP;
    result += newData.toHex();
    return result;
}

bool U2DbiPackUtils::unpackUdr(const QByteArray& modDetails, QByteArray& oldData, QByteArray& newData) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(tokens.count() == 2, QString("Invalid modDetails, wrong tokens count: %1. Expected - 2.").arg(tokens.size()), false);

    oldData = QByteArray::fromHex(tokens.first());
    newData = QByteArray::fromHex(tokens.last());

    return true;
}

QByteArray U2DbiPackUtils::packFeatureLocation(const U2DataId& featureId, const U2FeatureLocation &oldLocation, const U2FeatureLocation &newLocation) {
    QByteArray result;
    result += featureId.toHex();
    result += SEP;
    result += packOneFeatureLocation(oldLocation, SECOND_SEP);
    result += SEP;
    result += packOneFeatureLocation(newLocation, SECOND_SEP);
    return result;
}

bool U2DbiPackUtils::unpackFeatureLocation(const QByteArray &modDetails, U2DataId& featureId, U2FeatureLocation &oldLocation, U2FeatureLocation &newLocation) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    SAFE_POINT(tokens.count() == 3, QString("Invalid modDetails, wrong tokens count: %1. Expected - 3.").arg(tokens.size()), false);

    featureId = QByteArray::fromHex(tokens[0]);

    bool ok = unpackOneFeatureLocation(tokens[1], oldLocation, SECOND_SEP);
    CHECK(ok, false);

    ok = unpackOneFeatureLocation(tokens[2], newLocation, SECOND_SEP);
    CHECK(ok, false);

    return true;
}

QByteArray U2DbiPackUtils::packFeature(const U2Feature& oldFeature, const QList<U2FeatureKey>& oldKeys,
                              const U2Feature& newFeature, const QList<U2FeatureKey>& newKeys) {
    QByteArray result;
    result += packOneFeature(oldFeature, oldKeys);
    result += SEP;
    result += packOneFeature(newFeature, newKeys);
    return result;
}

bool U2DbiPackUtils::unpackFeature(const QByteArray& modDetails,
                          U2Feature& oldFeature, QList<U2FeatureKey>& oldKeys,
                          U2Feature& newFeature, QList<U2FeatureKey>& newKeys) {
    return true;
}

QByteArray U2DbiPackUtils::packOneFeature(const U2Feature& feature, const QList<U2FeatureKey>& keys) {
    QByteArray result; // 8 items

    result += feature.id.toHex();
    result += SEP;
    result += feature.sequenceId.toHex();
    result += SEP;
    result += feature.parentFeatureId.toHex();
    result += SEP;
    result += feature.rootFeatureId.toHex();
    result += SEP;
    result += feature.name;
    result += SEP;
    result += packOneFeatureLocation(feature.location, SECOND_SEP);
    result += SEP;
    result += QByteArray::number(feature.featureClass);
    result += SEP;
    result += QByteArray::number(feature.featureType);
    result += SEP;
    result += packFeatureKeys(keys); // pack keys number?

    return result;
}

bool U2DbiPackUtils::unpackOneFeature(const QByteArray& modDetails, U2Feature& feature, QList<U2FeatureKey>& keys) {
    QList<QByteArray> tokens = modDetails.split(SEP);
    //! TODO_SVEDIT: return back afters keys!
    SAFE_POINT(tokens.count() == 9, QString("Invalid modDetails, wrong tokens count: %1. Expected - 9.").arg(tokens.size()), false);

    feature.id = QByteArray::fromHex(tokens[0]);
    feature.sequenceId = QByteArray::fromHex(tokens[1]);
    feature.parentFeatureId = QByteArray::fromHex(tokens[2]);
    feature.rootFeatureId = QByteArray::fromHex(tokens[3]);
    feature.name = QString(tokens[4]);
    bool ok = unpackOneFeatureLocation(tokens[5], feature.location, SECOND_SEP);
    CHECK(ok, false);
    feature.featureClass = (U2Feature::FeatureClass)(tokens[6].toUInt(&ok)); //??
    CHECK(ok, false);
    feature.featureType = (U2FeatureType)(tokens[7].toUInt(&ok)); //??
    CHECK(ok, false);
    // TODO_SVEDIT: something is wrong here!
    ok = unpackFeatureKeys(tokens[8], keys);
    CHECK(ok, false);

    return true;
}

QByteArray U2DbiPackUtils::packOneFeatureLocation(const U2FeatureLocation& location, const char sep) {
    QByteArray result;
    result += QByteArray::number(location.strand.getDirectionValue());
    result += sep;
    result += QByteArray::number(location.region.startPos);
    result += sep;
    result += QByteArray::number(location.region.length);
    return result;
}

bool U2DbiPackUtils::unpackOneFeatureLocation(const QByteArray &modDetails, U2FeatureLocation& location, const char sep) {
    QList<QByteArray> tokens = modDetails.split(sep);
    SAFE_POINT(tokens.count() == 3, QString("Invalid modDetails, wrong tokens count: %1. Expected - 3.").arg(tokens.size()), false);

    bool ok = false;
    location.strand = U2Strand((U2Strand::Direction)tokens[0].toInt(&ok)); // TODO_SVEDIT: cast error
    CHECK(ok, false);
    qint64 start = tokens[1].toLongLong(&ok);
    CHECK(ok, false);
    qint64 len = tokens[2].toLongLong(&ok);
    CHECK(ok, false);
    location.region = U2Region(start, len);

    return true;
}

QByteArray U2DbiPackUtils::packFeatureKeys(const QList<U2FeatureKey>& keys) {
    QByteArray result;
    result += QByteArray::number(keys.size());
    result += SECOND_SEP;
    foreach (U2FeatureKey k, keys) {
        result += k.name;
        result += THIRD_SEP;
        result += k.value;
        result += SECOND_SEP;
    }
    return result;
}

bool U2DbiPackUtils::unpackFeatureKeys(const QByteArray &modDetails, QList<U2FeatureKey>& keys) {
    QList<QByteArray> tokens = modDetails.split(SECOND_SEP);
    int size = tokens[0].toInt();
    CHECK(size != 0, true);
    for (int i = 1; i < size + 1; i++) {
//    foreach (QByteArray t, tokens) {
        QByteArray t = tokens[i];
        QList<QByteArray> innerTokens = t.split(THIRD_SEP);
//        coreLog.info(QString(t.to));
        SAFE_POINT(innerTokens.count() == 2, QString("Invalid modDetails, wrong tokens count: %1. Expected - 2 for feature key.").arg(innerTokens.size()), false);
        U2FeatureKey featureKey(innerTokens.first(), innerTokens.last());
        keys.append(featureKey);
    }
    return true;
}

} // U2
