/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2016 UniPro <ugene@unipro.ru>
 * http://ugene.unipro.ru
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

#include <time.h>

#include <U2Core/DNAAlphabet.h>
#include <U2Core/DNAChromatogramObject.h>
#include <U2Core/DNAInfo.h>
#include <U2Core/DNASequenceObject.h>
#include <U2Core/GObjectReference.h>
#include <U2Core/GObjectRelationRoles.h>
#include <U2Core/GObjectTypes.h>
#include <U2Core/IOAdapter.h>
#include <U2Core/L10n.h>
#include <U2Core/TextObject.h>
#include <U2Core/TextUtils.h>
#include <U2Core/U2DbiUtils.h>
#include <U2Core/U2ObjectDbi.h>
#include <U2Core/U2OpStatus.h>
#include <U2Core/U2SafePoints.h>

#include "FsaFormat.h"
#include "DocumentFormatUtils.h"
#include "IOLibUtils.h"

/* TRANSLATOR U2::ABIFormat */

namespace U2 {

FsaFormat::FsaFormat(QObject* p) : DocumentFormat(p, DocumentFormatFlag_SupportStreaming, QStringList() << "fsa") {
    formatName = tr("FSA");
    formatDescription = tr("A chromatogram file format");
    supportedObjectTypes += GObjectTypes::TEXT;
}

FormatCheckResult FsaFormat::checkRawData(const QByteArray& rawData, const GUrl&) const {
    const char* data = rawData.constData();
    int size = rawData.size();

    if (size <= 4 || data[0] != 'A' || data[1] != 'B' || data[2] != 'I' || data[3] != 'F') {
        /*
        * Maybe we've got a file in MacBinary format in which case
        * we'll have an extra offset 128 bytes to add.
        */
        data += 128;
        size -= 128;
        if (size <= 4 || data[0] != 'A' || data[1] != 'B' || data[2] != 'I' || data[3] != 'F') {
            return FormatDetection_NotMatched;
        }
    }
    bool hasBinaryBlocks = TextUtils::contains(TextUtils::BINARY, data, size);
    return hasBinaryBlocks ? FormatDetection_Matched : FormatDetection_NotMatched;
}

Document* FsaFormat::loadDocument(IOAdapter* io, const U2DbiRef& dbiRef, const QVariantMap& fs, U2OpStatus& os) {
    QByteArray readBuff;
    QByteArray block(BUFF_SIZE, 0);
    quint64 len = 0;
    while ((len = io->readBlock(block.data(), BUFF_SIZE)) > 0) {
        readBuff.append(QByteArray(block.data(), len));
        CHECK_EXT(readBuff.size() <= CHECK_MB, os.setError(L10N::errorFileTooLarge(io->getURL())), NULL);
    }

    SeekableBuf sf;
    sf.head = readBuff.constData();
    sf.pos = 0;
    sf.size = readBuff.size();
    Document* doc = parseABI(dbiRef, &sf, io, fs, os);
    CHECK_OP(os, NULL)
        CHECK_EXT(doc != NULL, os.setError(tr("Not a valid ABIF file: %1").arg(io->toString())), NULL);
    return doc;
}

/*
* Copyright (c) Medical Research Council 1994. All rights reserved.
*
* Permission to use, copy, modify and distribute this software and its
* documentation for any purpose is hereby granted without fee, provided that
* this copyright and notice appears in all copies.
*
* This file was written by James Bonfield, Simon Dear, Rodger Staden,
* as part of the Staden Package at the MRC Laboratory of Molecular
* Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
*
* MRC disclaims all warranties with regard to this software.
*/

/*
* Title:    seqIOABI
*
* File:     seqIOABI.c
* Purpose:  Reading (not writing) of ABI sequences
* Last update: Fri Sep 02, 1994
*/

/*
* The ABI magic number - "ABIF"
*/
#define ABI_MAGIC    ((int) ((((('A'<<8)+'B')<<8)+'I')<<8)+'F')

/*
* The index is located towards the end of the ABI trace file.
* It's location is given by a longword at a fixed place.
*/
#define IndexPO 26

#define IndexEntryLength 28

/*
* Here are some labels we will be looking for, four chars packed
* into an int_4
*/
#define LABEL(a) ((int) ((((((a)[0]<<8)+(a)[1])<<8)+(a)[2])<<8)+(a)[3])

//pStrings
#define CMNTLabel         LABEL("CMNT")
#define LIMSLabel         LABEL("LIMS")
#define SpNmLabel         LABEL("SpNm")
#define StdFLabel         LABEL("StdF")
#define DySNLabel         LABEL("DySN")
#define GTypLabel         LABEL("GTyp")
#define MCHNLabel         LABEL("MCHN")
#define SMLtLabel         LABEL("SMLt")
#define SVERLabel         LABEL("SVER")
#define TUBELabel         LABEL("TUBE")
#define USERLabel         LABEL("USER")
#define DyeNLabel         LABEL("DyeN")
#define HCFGLabel         LABEL("HCFG")

//cStrings
#define STYPLabel         LABEL("STYP")
#define CTIDLabel         LABEL("CTID")
#define RunNLabel         LABEL("RunN")

//long
#define EPVtLabel         LABEL("EPVt")
#define InScLabel         LABEL("InSc")
#define InVtLabel         LABEL("InVt")
#define LsrPLabel         LABEL("LsrP")
#define PSZELabel         LABEL("PSZE")
#define PXLBLabel         LABEL("PXLB")
#define TmprLabel         LABEL("Tmpr")

//short
#define DSamLabel         LABEL("DSam")
#define LANELabel         LABEL("LANE")
#define LNTDLabel         LABEL("LNTD")
#define NLNELabel         LABEL("NLNE")
#define DyeWLabel         LABEL("DyeW")

//char
#define CpEPLabel         LABEL("CpEP")
#define MODLLabel         LABEL("MODL")

//date+time
#define RUNDLabel         LABEL("RUNDLabel")
#define RUNTLabel         LABEL("RUNTLabel")

//data
#define DataEntryLabel    LABEL("DATA")


/*
#define DataEntryLabel    LABEL("DATA")
#define FWO_Label         LABEL("FWO_")

#define PDMFLabel         LABEL("PDMF")
#define SMPLLabel         LABEL("SMPL")

#define LANELabel         LABEL("LANE")
#define RUNDLabel         LABEL("RUND")
#define RUNTLabel         LABEL("RUNT")
#define SPACLabel         LABEL("SPAC")
#define SVERLabel         LABEL("SVER")
#define MODLLabel         LABEL("MODL")
*/

#define baseIndex(B) ((B)=='C'?0:(B)=='A'?1:(B)=='G'?2:3)
namespace {
/*
* Gets the offset of the ABI index.
* Returns -1 for failure, 0 for success.
*/
    int getABIIndexOffset(SeekableBuf* fp, uint *indexO) {
        uint magic = 0;

        /*
        * Initialise header_fudge.
        *
        * This is usually zero, but maybe we've transfered a file in MacBinary
        * format in which case we'll have an extra 128 bytes to add to all
        * our fseeks.
        */
        be_read_int_4(fp, &magic);
        if (magic != ABI_MAGIC) {
            fp->head += 128;
            fp->size -= 128;
        }

        if ((SeekBuf(fp, IndexPO, 0) != 0) || (!be_read_int_4(fp, indexO)))
            return -1;
        else
            return 0;
    }

    /*
    * From the ABI results file connected to `fp' whose index starts
    * at byte offset `indexO', return in `val' the `lw'th long word
    * from the `count'th entry labelled `label'.
    * The result is 0 for failure, or index offset for success.
    */
    int getABIIndexEntryLW(SeekableBuf* fp, int indexO, uint label, uint count, int lw, uint *val) {
        int entryNum = -1;
        int i;
        uint entryLabel, entryLw1;

        do {
            entryNum++;

            if (SeekBuf(fp, indexO + (entryNum*IndexEntryLength), 0) != 0)
                return 0;

            if (!be_read_int_4(fp, &entryLabel))
                return 0;

            if (!be_read_int_4(fp, &entryLw1))
                return 0;
        } while (!(entryLabel == label && entryLw1 == count));

        for (i = 2; i <= lw; i++) {
            if (!be_read_int_4(fp, val))
                return 0;
        }

        return indexO + (entryNum*IndexEntryLength);
    }

    int cStringEndPos(SeekableBuf* fp, int indexO, uint label) {
        int entryNum = -1;
        ushort entryLw1;

        while (true) {
            entryNum++;
            if (SeekBuf(fp, indexO + (entryNum), 0) != 0)
                return 0;
            be_read_int_2(fp, &entryLw1);
            if (entryLw1 == 0)
                return entryNum;
        }
        return 0;
    }

    /*
    * From the ABI results file connected to `fp' whose index starts
    * at byte offset `indexO', return in `val' the `sw'th short word
    * from the `count'th entry labelled `label'.
    * The result is 0 for failure, or index offset for success.
    */
    int getABIIndexEntrySW(SeekableBuf* fp, int indexO, uint label, uint count, int sw, ushort *val) {
        int entryNum = -1;
        int i;
        uint entryLabel, entryLw1, position;

        do {
            entryNum++;
            position = indexO + (entryNum*IndexEntryLength);
            if (SeekBuf(fp, position, 0) != 0)
                return 0;

            if (!be_read_int_4(fp, &entryLabel))
                return 0;

            if (!be_read_int_4(fp, &entryLw1))
                return 0;
        } while (!(entryLabel == label && entryLw1 == count));

        for (i = 4; i <= sw; i++) {
            if (!be_read_int_2(fp, val))
                return 0;
        }

        return indexO + (entryNum*IndexEntryLength);
    }

    /*
    * Get an "ABI String". These strings are either pointed to by the index
    * offset, or held in the offset itself when the string is <= 4 characters.
    * The "type" of the index entry is either 0x12 (a pascal string in which
    * case the first byte of the string determines its length) or a 0x02 (a
    * C-style string with length coming from the abi index).
    *
    * "string" will be max 256 bytes for the pascal string, but is of unknown
    * (and hence potentially buggy) length for C-strings. For now we live with
    * it as this entire file needs rewriting from scratch anyway.
    *
    * Returns -1 for failure, string length for success.
    */
    int getABIString(SeekableBuf *fp, int indexO, uint label, uint count, char *string) {
        uint off;
        uint len;
        quint16 type;

        off = getABIIndexEntrySW(fp, indexO, label, count, 4, &type);
        if (!off)
            return -1;

        if ((off = getABIIndexEntryLW(fp, indexO, label, count, 4, &len))) {
            uchar len2 = 0;

            if (!len)
                return 0;

            /* Determine offset */
            if (len <= 4)
                off += 20;
            else
                getABIIndexEntryLW(fp, indexO, label, count, 5, &off);

            /* Read length byte */
            if (type == 0x12) {
                SeekBuf(fp, off, 0);
                be_read_int_1(fp, &len2);
            } else {
                len2 = len;
            }

            /* Read data (max 255 bytes) */
            fp->read(string, len2);
            string[len2] = 0;

            return len2;
        } else {
            return -1;
        }
    }

    int getCString(SeekableBuf *fp, int indexO, uint label, uint count, char *string) {
        uint off;
        uint len;
        quint16 type;

        off = getABIIndexEntrySW(fp, indexO, label, count, 4, &type);
        if (!off)
            return -1;

        if ((off = getABIIndexEntryLW(fp, indexO, label, count, 4, &len))) {
            if (!len)
                return 0;

            /* Determine offset */
            /*
            if (len <= 4)
                off += 20;
            else
            */
            getABIIndexEntryLW(fp, indexO, label, count, 5, &off);
            
            SeekBuf(fp, off, 0);
            /* Read data (max 255 bytes) */
            fp->read(string, len);
            string[len] = 0;

            return len;
        } else {
            return -1;
        }
    }

    /*
    * Get an "ABI Int_1". This is raw 1-byte integer data pointed to by the
    * offset, or held in the offset itself when the data is <= 4 characters.
    *
    * If indexO is 0 then we do not search for (or indeed use) label and count,
    * but simply assume that we are already at the correct offset and read from
    * here. (NB: This negates the length <= 4 check.)
    *
    * Returns -1 for failure, length desired for success (it'll only fill out
    * up to max_data_len elements, but it gives an indication of whether there
    * was more to come).
    */
    int getABIint1(SeekableBuf *fp, int indexO, uint label, uint count, uchar *data, int max_data_len) {
        uint off;
        uint len, len2;

        if (indexO) {
            if (!(off = getABIIndexEntryLW(fp, indexO, label, count, 4, &len)))
                return -1;

            if (!len)
                return 0;

            /* Determine offset */
            if (len <= 4)
                off += 20;
            else
                getABIIndexEntryLW(fp, indexO, label, count, 5, &off);

            len2 = qMin((uint)max_data_len, len);

            SeekBuf(fp, off, 0);
        } else {
            len = len2 = max_data_len;
        }

        fp->read((char*)data, len2);

        return len;
    }

    /*
    * Get an "ABI Int_2". This is raw 2-byte integer data pointed to by the
    * offset, or held in the offset itself when the data is <= 4 characters.
    *
    * Returns -1 for failure, length desired for success (it'll only fill out
    * up to max_data_len elements, but it gives an indication of whether there
    * was more to come).
    */
    int getABIint2(SeekableBuf *fp, int indexO, uint label, uint count, ushort *data, int max_data_len) {
        int len, l2;
        int i;

        len = getABIint1(fp, indexO, label, count, (uchar *)data, max_data_len * 2);
        if (-1 == len)
            return -1;

        len /= 2;
        l2 = qMin(len, max_data_len);
        for (i = 0; i < l2; i++) {
            data[i] = be_int2((uchar*)(data + i));
        }

        return len;
    }

    /*
    * Get an "ABI Int_4". This is raw 4-byte integer data pointed to by the
    * offset, or held in the offset itself when the data is <= 4 characters.
    *
    * Returns -1 for failure, length desired for success (it'll only fill out
    * up to max_data_len elements, but it gives an indication of whether there
    * was more to come).
    */
    int getABIint4(SeekableBuf *fp, int indexO, uint label, uint count, uint *data, int max_data_len) {
        int len, l2;
        int i;

        len = getABIint1(fp, indexO, label, count, (uchar *)data, max_data_len * 4);
        if (-1 == len)
            return -1;

        len /= 4;
        l2 = qMin(len, max_data_len);
        for (i = 0; i < l2; i++) {
            data[i] = be_int4((uchar*)(data + i));
        }

        return len;
    }

    void replace_nl(char *string) {
        char *cp;

        for (cp = string; *cp; cp++) {
            if (*cp == '\n') *cp = ' ';
        }
    }
}

Document* FsaFormat::parseABI(const U2DbiRef& dbiRef, SeekableBuf* fp, IOAdapter* io, const QVariantMap& fs, U2OpStatus& os) {
    DbiOperationsBlock opBlock(dbiRef, os);
    CHECK_OP(os, NULL);
    Q_UNUSED(opBlock);
    QList<GObject*> objects;
    QString textData;

    if (!loadTextObject(fp, textData)) {
        return NULL;
    }

    QVariantMap hints;
    hints.insert(DBI_FOLDER_HINT, fs.value(DBI_FOLDER_HINT, U2ObjectDbi::ROOT_FOLDER));
    TextObject *textObj = TextObject::createInstance(textData, "Info", dbiRef, os, hints);
    CHECK_OP(os, NULL);
    objects.append(textObj);

    Document* doc = new Document(this, io->getFactory(), io->getURL(), dbiRef, objects, fs);
    return doc;
}

bool FsaFormat::loadTextObject(SeekableBuf* fp, QString &textData) {
    uint indexO;

    if (-1 == getABIIndexOffset(fp, &indexO)) {
        return false;
    }
    //pStrings
    extactPString(fp, indexO, CMNTLabel, textData, QString("Comment: %1\n"));
    extactPString(fp, indexO, LIMSLabel, textData, QString("Tracking ID: %1\n"));
    extactPString(fp, indexO, SpNmLabel, textData, QString("Name: %1\n"));
    extactPString(fp, indexO, StdFLabel, textData, QString("Size Standard name: %1\n"));
    extactPString(fp, indexO, DySNLabel, textData, QString("Dye set name: %1\n"));
    extactPString(fp, indexO, GTypLabel, textData, QString("Gel type: %1\n"));
    extactPString(fp, indexO, MCHNLabel, textData, QString("Machine name: %1\n"));    
    extactPString(fp, indexO, SpNmLabel, textData, QString("Separation member lot number: %1\n"));
    extactPString(fp, indexO, SVERLabel, textData, QString("Data collection version number: %1\n"));
    extactPString(fp, indexO, SVERLabel, textData, QString("Firmware version number: %1\n"), 3);
    extactPString(fp, indexO, TUBELabel, textData, QString("Autosampler position: %1\n"));
    extactPString(fp, indexO, USERLabel, textData, QString("User name of plate creator: %1\n"));
    extactPString(fp, indexO, DyeNLabel, textData, QString("Dye names: %1\n"));
    
    //cStrings
    extactCString(fp, indexO, STYPLabel, textData, QString("Type: %1\n"));
    extactCString(fp, indexO, CTIDLabel, textData, QString("Container ID: %1\n"));
    extactCString(fp, indexO, RunNLabel, textData, QString("Run name: %1\n"));
    extactCString(fp, indexO, HCFGLabel, textData, QString("Instrument parameters: %1\n"), 4);
    
    //long
    extractLong(fp, indexO, EPVtLabel, textData, QString("Electrophoresis voltage: %1\n"));
    extractLong(fp, indexO, InScLabel, textData, QString("Injection time: %1\n"));
    extractLong(fp, indexO, InVtLabel, textData, QString("Injection voltage: %1\n"));
    extractLong(fp, indexO, LsrPLabel, textData, QString("Laser power setting in micro Watts: %1\n"));
    extractLong(fp, indexO, PSZELabel, textData, QString("Plate size: %1\n"));
    extractLong(fp, indexO, PXLBLabel, textData, QString("Pixel bin size: %1\n"));
    extractLong(fp, indexO, TmprLabel, textData, QString("Oven temperature: %1\n"));

    //short
    extractShort(fp, indexO, DSamLabel, textData, QString("Downsampling rate: %1\n"));
    extractShort(fp, indexO, LANELabel, textData, QString("Sample's lane or capillary number: %1\n"));
    extractShort(fp, indexO, LNTDLabel, textData, QString("Length to detector: %1\n"));
    extractShort(fp, indexO, NLNELabel, textData, QString("Total number of capillaries : %1\n"));
    extractShort(fp, indexO, DyeWLabel, textData, QString("Dye wawelengths: %1\n"));

    //char
    extractShort(fp, indexO, CpEPLabel, textData, QString("Capillary type electrophoresis: %1\n"));
    extractShort(fp, indexO, MODLLabel, textData, QString("Model number: %1\n"));

    //date-time
    extractDateTime(fp, indexO, textData);

    //read data arrays
    extractData(fp, indexO, textData);

    return true;
}

void FsaFormat::extactPString(SeekableBuf *fp, uint indexO, uint label, QString &textData, const QString &stringPattern, uint entryIndex /*= 1*/) {
    char strBuf[256];
    int clen = getABIString(fp, indexO, label, entryIndex, strBuf);
    if (clen != -1) {
        strBuf[clen] = 0;
        char *commstrp = strBuf;
        char *p;
        do {
            if ((p = strchr(commstrp, '\n'))) {
                *p++ = 0;
            }
            textData.append(stringPattern.arg(commstrp));
        } while ((commstrp = p));
    }
}

void FsaFormat::extactCString(SeekableBuf *fp, uint indexO, uint label, QString &textData, const QString &stringPattern, uint entryIndex /*= 1*/) {
    char strBuf[256];
    int clen = getCString(fp, indexO, label, entryIndex, strBuf);
    if (clen != -1) {
        strBuf[clen] = 0;
        char *commstrp = strBuf;
        char *p;
        do {
            if ((p = strchr(commstrp, '\n'))) {
                *p++ = 0;
            }
            textData.append(stringPattern.arg(commstrp));
        } while ((commstrp = p));
    }
}

void FsaFormat::extractLong(SeekableBuf* fp, uint indexO, uint label, QString &textData, const QString &stringPattern) {
    ulong longVal = 0;
    if (-1 != getABIint4(fp, indexO, label, 1, (uint *)&longVal, 1)) {
        textData.append(stringPattern.arg(longVal));
    }
}

void FsaFormat::extractShort(SeekableBuf* fp, uint indexO, uint label, QString &textData, const QString &stringPattern) {
    ushort shortVal = 0;
    if (-1 != getABIint2(fp, indexO, label, 1, (ushort *)&shortVal, 1)) {
        textData.append(stringPattern.arg(shortVal));
    }
}

void FsaFormat::extractDateTime(SeekableBuf* fp, uint indexO, QString &textData) {
    uint offset = 0;
    uint offset2 = 0;
    uint offset3 = 0;
    uint offset4 = 0;
    if (getABIIndexEntryLW(fp, indexO, RUNDLabel, 1, 5, &offset) &&
        getABIIndexEntryLW(fp, indexO, RUNDLabel, 2, 5, &offset2) &&
        getABIIndexEntryLW(fp, indexO, RUNTLabel, 1, 5, &offset3) &&
        getABIIndexEntryLW(fp, indexO, RUNTLabel, 2, 5, &offset4)) {
        //char buffer[1025];
        char buffer_s[1025];
        char buffer_e[1025];
        struct tm t;
        uint rund_s, rund_e, runt_s, runt_e;

        rund_s = offset;
        rund_e = offset2;
        runt_s = offset3;
        runt_e = offset4;

        QString buffer = QString("%1%2%3.%4%5%6 - %7%8%9.%10%11%12")
            .arg((rund_s >> 16), 4, 10, QLatin1Char('0')).arg((rund_s >> 8) & 0xff, 2, 10, QLatin1Char('0'))
            .arg((rund_s & 0xff), 2, 10, QLatin1Char('0')).arg(runt_s >> 24, 2, 10, QLatin1Char('0'))
            .arg((runt_s >> 16) & 0xff, 2, 10, QLatin1Char('0')).arg((runt_s >> 8) & 0xff, 2, 10, QLatin1Char('0'))
            .arg(rund_e >> 16, 4, 10, QLatin1Char('0')).arg((rund_e >> 8) & 0xff, 2, 10, QLatin1Char('0'))
            .arg(rund_e & 0xff, 2, 10, QLatin1Char('0')).arg(runt_e >> 24, 2, 10, QLatin1Char('0'))
            .arg((runt_e >> 16) & 0xff, 2, 10, QLatin1Char('0')).arg((runt_e >> 8) & 0xff, 2, 10, QLatin1Char('0'));

        memset(&t, 0, sizeof(t));
        t.tm_mday = rund_s & 0xff;
        t.tm_mon = ((rund_s >> 8) & 0xff) - 1;
        t.tm_year = (rund_s >> 16) - 1900;
        t.tm_hour = runt_s >> 24;
        t.tm_min = (runt_s >> 16) & 0xff;
        t.tm_sec = (runt_s >> 8) & 0xff;
        t.tm_isdst = -1;
        /*
        * Convert struct tm to time_t. We ignore the time_t value, but
        * the conversion process will update the tm_wday element of
        * struct tm.
        */
        mktime(&t);
        strftime(buffer_s, 1024, "%a %d %b %H:%M:%S %Y", &t);

        t.tm_mday = rund_e & 0xff;
        t.tm_mon = ((rund_e >> 8) & 0xff) - 1;
        t.tm_year = (rund_e >> 16) - 1900;
        t.tm_hour = runt_e >> 24;
        t.tm_min = (runt_e >> 16) & 0xff;
        t.tm_sec = (runt_e >> 8) & 0xff;
        t.tm_isdst = -1;
        /*
        * Convert struct tm to time_t. We ignore the time_t value, but
        * the conversion process will update the tm_wday element of
        * struct tm.
        */
        mktime(&t);
        strftime(buffer_e, 1024, "%a %d %b %H:%M:%S %Y", &t);

        textData.append(QString("DATE=%1 to %2\nRUND=%3\n").arg(buffer_s).arg(buffer_e).arg(buffer));
    }
}

void FsaFormat::extractData(SeekableBuf* fp, uint indexO, QString &textData) {
    uint firstElementOffset, arraySize;
    int result1, result2;
    textData.append("================= CSV data by points below =================\n");
    for (int i = 1; i <= 4; i++) {
        result1 = getABIIndexEntryLW(fp, indexO, DataEntryLabel, i, 5, &firstElementOffset);
        result2 = getABIIndexEntryLW(fp, indexO, DataEntryLabel, i, 3, &arraySize);
        if (result1 != 0 && result2 != 0) {
            QVector<int> vec;
            SeekBuf(fp, firstElementOffset, 0);
            for (uint itemNumber = 0; itemNumber < arraySize; itemNumber ++) {
                char strBuf[2];
                fp->read(strBuf, 2);
                short readedData = ((strBuf[0] << 8) & 0xFF00) | (strBuf[1] & 0xFF);
                vec.append(readedData);
            }
            dataMap.insert(i, vec);
        }
    }

    for (int i = 105; i <= 199; i++) {
        result1 = getABIIndexEntryLW(fp, indexO, DataEntryLabel, i, 5, &firstElementOffset);
        result2 = getABIIndexEntryLW(fp, indexO, DataEntryLabel, i, 3, &arraySize);
        if (result1 != 0 && result2 != 0) {
            QVector<int> vec;
            SeekBuf(fp, firstElementOffset, 0);
            for (uint itemNumber = 0; itemNumber < arraySize; itemNumber++) {
                char strBuf[2];
                fp->read(strBuf, 2);
                short readedData = ((strBuf[0] << 8) & 0xFF00) | (strBuf[1] & 0xFF);
                vec.append(readedData);
            }
            dataMap.insert(i, vec);
        }
    }

    foreach(int key, dataMap.keys()) {
        textData.append(QString::number(key));
        textData.append(";");
    }
    textData.append("\n");
    if (dataMap.isEmpty()) {
        return;
    }
    for (int i = 0; i < dataMap.first().size(); i++) {
        foreach(int key, dataMap.keys()) {
            textData.append(QString::number(dataMap[key][i]));
            textData.append(";");
        }
        textData.append("\n");
    }
}

}//namespace
