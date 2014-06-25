/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2014 UniPro <ugene@unipro.ru>
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

#include <QtCore/QScopedArrayPointer>
#include <QtCore/QScopedPointer>

#include <U2Core/AppContext.h>
#include <U2Core/AnnotationTableObject.h>
#include <U2Core/GObjectReference.h>
#include <U2Core/GObjectRelationRoles.h>
#include <U2Core/IOAdapter.h>
#include <U2Core/Log.h>
#include <U2Core/L10n.h>
#include <U2Core/TextUtils.h>
#include <U2Core/U2SafePoints.h>
#include <U2Core/U1AnnotationUtils.h>
#include <U2Core/U2DbiRegistry.h>
#include <U2Core/U2DbiUtils.h>
#include <U2Core/U2ObjectDbi.h>
#include <U2Core/U2OpStatus.h>
#include <U2Core/U2SafePoints.h>
#include <U2Core/U2SequenceUtils.h>

#include "DocumentFormatUtils.h"
#include "BedFormat.h"

namespace U2{

//-------------------------------------------------------------------
//  BEDLineValidateFlags
//-------------------------------------------------------------------
BEDLineValidateFlags::BEDLineValidateFlags()
    : incorrectNumberOfFields(false),
      emptyFields(false),
      incorrectCoordinates(false),
      incorrectScore(false),
      incorrectStrand(false),
      incorrectThickCoordinates(false),
      incorrectItemRgb(false),
      incorrectBlocks(false),
      hasTrackLine(false)
{
}

FormatDetectionScore BEDLineValidateFlags::getFormatDetectionScore()
{
    if (incorrectNumberOfFields || emptyFields || incorrectCoordinates) {
        return FormatDetection_NotMatched;
    }

    if (incorrectScore || incorrectStrand || incorrectBlocks || !hasTrackLine) {
        return FormatDetection_LowSimilarity;
    }

    if (incorrectThickCoordinates || incorrectItemRgb) {
        return FormatDetection_HighSimilarity;
    }

    return FormatDetection_Matched;
}


//-------------------------------------------------------------------
//  BedFormat
//-------------------------------------------------------------------
const QString BedFormat::FORMAT_NAME = BedFormat::tr("BED");

const QString BedFormat::TRACK_NAME_QUALIFIER_NAME = "track_name";
const QString BedFormat::TRACK_DESCR_QUALIFIER_NAME = "track_description";
const QString BedFormat::CHROM_QUALIFIER_NAME = "chrom";
const QString BedFormat::ANNOT_QUALIFIER_NAME = "name";
const QString BedFormat::SCORE_QUALIFIER_NAME = "score";
const QString BedFormat::STRAND_QUALIFIER_NAME = "strand";
const QString BedFormat::THICK_START_QUALIFIER_NAME = "thick_start";
const QString BedFormat::THICK_END_QUALIFIER_NAME = "thick_end";
const QString BedFormat::ITEM_RGB_QUALIFIER_NAME = "item_rgb";
const QString BedFormat::BLOCK_COUNT_QUALIFIER_NAME = "block_count";
const QString BedFormat::BLOCK_SIZES_QULAIFIER_NAME = "block_sizes";
const QString BedFormat::BLOCK_STARTS_QUALIFIER_NAME = "block_starts";


BedFormat::BedFormat(QObject* p)
    : DocumentFormat(p, DocumentFormatFlag_SupportWriting, QStringList("bed"))
{
    formatDescription = tr("The BED (Browser Extensible Data) format was developed by UCSC for displaying transcript structures in the genome browser.");
    supportedObjectTypes+=GObjectTypes::ANNOTATION_TABLE;
}


Document* BedFormat::loadDocument(IOAdapter* io, const U2DbiRef& dbiRef, const QVariantMap& fs, U2OpStatus& os)
{
    CHECK_EXT(io != NULL && io->isOpen(), os.setError(L10N::badArgument("IO adapter")), NULL);
    QList<GObject*> objects;

    load(io, objects, dbiRef, os, fs);
    CHECK_OP_EXT(os, qDeleteAll(objects), NULL);

    Document* doc = new Document(this, io->getFactory(), io->getURL(), dbiRef, objects);
    return doc;
}

void BedFormat::load(IOAdapter* io, QList<GObject*>& objects, const U2DbiRef& dbiRef, U2OpStatus& os, const QVariantMap& fs) {
    QString defaultAnnotName = "misc_feature";
    const QHash<QString, QList<SharedAnnotationData> >& annotationsHash = parseDocument(io, defaultAnnotName, os);

    foreach ( const QString &sequenceName, annotationsHash.keys( ) ) {
        const QString annotTableName = sequenceName + FEATURES_TAG;
        AnnotationTableObject *annotTable = NULL;
        foreach (GObject* object, objects) {
            if (object->getGObjectName() == annotTableName) {
                annotTable = dynamic_cast<AnnotationTableObject *>( object );
            }
        }
        if (!annotTable) {
            QVariantMap hints;
            hints.insert(DBI_FOLDER_HINT, fs.value(DBI_FOLDER_HINT, U2ObjectDbi::ROOT_FOLDER));
            annotTable = new AnnotationTableObject( annotTableName, dbiRef, hints);
            objects.append(annotTable);
        }

        // Assume that the group name is the same as the annotation name
        QString groupName = defaultAnnotName;
        if (AnnotationGroup::isValidGroupName(groupName, false)) {
            groupName = "Group"; // or set this name if the annotation name is not appropriate
        }

        const QList<SharedAnnotationData>& annotations = annotationsHash.value(sequenceName);
        foreach (const SharedAnnotationData& annotData, annotations) {
            annotTable->addAnnotation( *annotData , groupName);
        }
    }
}

/** Validate the values: they must be integer and within the bound of the region */
bool validateThickCoordinates(const QString& thickStartStr, const QString& thickEndStr, const U2Region& region)
{
    if (thickStartStr.isEmpty() || thickEndStr.isEmpty()) {
        return false;
    }

    bool thickStartIsInt;
    bool thickEndIsInt;
    int thickStart = thickStartStr.toInt(&thickStartIsInt);
    int thickEnd = thickEndStr.toInt(&thickEndIsInt);

    if (!thickStartIsInt || !thickEndIsInt ||
        (thickStart < region.startPos) ||
        (thickEnd > region.endPos()) ||
        (thickStart >= thickEnd)) {
            return false;
    }

    return true;
}


/**
* Validate that the color is valid
* If color is "0", then it is not set, otherwise it has format r,g,b
* If color is valid, corresponding QColor is set.
*/
bool validateAnnotationColor(const QString& itemRgbStr, QColor& annotColor)
{
    if ("0" == itemRgbStr) {
        return true;
    }

    QStringList rgbValues = itemRgbStr.split(",");
    if (3 != rgbValues.count()) {
        return false;
    }

    bool convertionStatusIsOk;

    int red = rgbValues[0].toInt(&convertionStatusIsOk);
    if (!convertionStatusIsOk) {
        return false;
    }

    int green = rgbValues[1].toInt(&convertionStatusIsOk);
    if (!convertionStatusIsOk) {
        return false;
    }

    int blue = rgbValues[2].toInt(&convertionStatusIsOk);
    if (!convertionStatusIsOk) {
        return false;
    }

    QColor color(red, green, blue);
    if (!color.isValid()) {
        return false;
    }

    annotColor = color;
    return true;
}


/**
 * Validate blocks (exons):
 * blockCount is the number of blocks.
 * blockSizes is a comma-separated list of the block sizes. The number of items must correspond to blockCount.
 * blockStarts is a comma-separated list of the block starts. All of the blockStart positions are calculated
 * relative to chromStart. The number of items must correspond to blockCount.
 * Note that the validated/calculated exons are not currently remembered.
 */
bool validateBlocks(const QString& blockCountStr, const QString& blockSizesStr, const QString& blockStartsStr, const U2Region& region)
{
    bool blockCountIsInt;
    int blockCount = blockCountStr.toInt(&blockCountIsInt);
    if (!blockCountIsInt || (blockCount == 0)) {
        return false;
    }

    // Skipping empty parts because values can have e.g. the following format: "567, 488," (comma at the end)
    QStringList blockSizesStrValues = blockSizesStr.split(",", QString::SkipEmptyParts);
    if (blockSizesStrValues.count() != blockCount) {
        return false;
    }

    QStringList blockStartsStrValues = blockStartsStr.split(",", QString::SkipEmptyParts);
    if (blockStartsStrValues.count() != blockCount) {
        return false;
    }

    bool conversionIsOk;
    QVector<int> blockSizes;
    QVector<int> blockStarts;
    for (int i = 0; i < blockCount; ++i) {
        blockSizesStrValues[i].toInt(&conversionIsOk);
        if (false == conversionIsOk) {
            return false;
        }

        int start = blockStartsStrValues[i].toInt(&conversionIsOk);
        if (false == conversionIsOk) {
            return false;
        }

        if (start > region.length) {
            return false;
        }
    }

    return true;
}


#define CHECK_FIELD(INDEX) \
    if (fields.size()<INDEX || INDEX<0) { \
        status.incorrectNumberOfFields = true; \
        return parsedData; \
    }

BedLineData BedFormat::parseAndValidateLine(const QString& line, int numOfFields, BEDLineValidateFlags& status) const
{
    BedLineData parsedData;

    // All fields are separated by a single tab
    QStringList fields = line.split("\t");

    // Verify that the line has the expected number of fields and they are not empty
    // Do not continue to validate the line if it is incorrect
    if (numOfFields != fields.count()) {
        status.incorrectNumberOfFields = true;
        return parsedData;
    }

    foreach (QString field, fields) {
        if (field.trimmed().isEmpty()) {
            status.emptyFields = true;
            return parsedData;
        }
    }

    // Coordinates
    // "start" can be zero, "end" is not included into the region
    bool startIsInt;
    bool endIsInt;
    CHECK_FIELD(BED_CHROM_START_INDEX);
    CHECK_FIELD(BED_CHROM_END_INDEX);
    qint64 start = fields[BED_CHROM_START_INDEX].toLongLong(&startIsInt);
    qint64 end = fields[BED_CHROM_END_INDEX].toLongLong(&endIsInt);
    if (!startIsInt || !endIsInt || (start < 0) || (start >= end)) {
        status.incorrectCoordinates = true;
        return parsedData;
    }

    // Fill in the data and continue validation even if a value is incorrect
    CHECK_FIELD(BED_CHROM_NAME_INDEX);
    parsedData.seqName = fields[BED_CHROM_NAME_INDEX];
    parsedData.region = U2Region(start, end - start);

    // Annotation name
    if (numOfFields > BED_ANNOT_NAME_INDEX) {
        CHECK_FIELD(BED_ANNOT_NAME_INDEX);
        parsedData.additionalFields[ANNOT_QUALIFIER_NAME] = fields[BED_ANNOT_NAME_INDEX];
    }

    // Score
    if (numOfFields > BED_SCORE_INDEX) {
        CHECK_FIELD(BED_SCORE_INDEX);
        QString scoreStr = fields[BED_SCORE_INDEX];
        parsedData.additionalFields[SCORE_QUALIFIER_NAME] = scoreStr;

        // Validate the value: it should be an integer value between 0 and 1000
        // UPDATE: for MACS peaks its wrong, so validation is removed
        bool scoreIsOk;
        scoreStr.toInt(&scoreIsOk);
        if (!scoreIsOk) {
            scoreStr.toDouble(&scoreIsOk);
        }
        if (!scoreIsOk) {
            status.incorrectScore = true;
        }
    }

    // Strand (either '+' or '-')
    if (numOfFields > BED_STRAND_INDEX) {
        CHECK_FIELD(BED_STRAND_INDEX);
        QString strandStr = fields[BED_STRAND_INDEX];
        parsedData.additionalFields[STRAND_QUALIFIER_NAME] = strandStr;

        if ("+" != strandStr &&
            "-" != strandStr) {
                status.incorrectStrand = true;
        }
    }

    // Thick coordinates
    if (numOfFields > BED_THICK_START_INDEX) {
        CHECK_FIELD(BED_THICK_START_INDEX);
        QString thickStartStr = fields[BED_THICK_START_INDEX];
        QString thickEndStr;

        // If thick start is set, thick end must also be set
        if (numOfFields <= BED_THICK_END_INDEX) {
            status.incorrectThickCoordinates = true;
        }
        else {
            CHECK_FIELD(BED_THICK_END_INDEX);
            thickEndStr = fields[BED_THICK_END_INDEX];
        }

        parsedData.additionalFields[THICK_START_QUALIFIER_NAME] = thickStartStr;
        parsedData.additionalFields[THICK_END_QUALIFIER_NAME] = thickEndStr;
    }

    // Annotation color
    if (numOfFields > BED_ITEM_RGB_INDEX) {
        CHECK_FIELD(BED_ITEM_RGB_INDEX);
        QString itemRgbStr = fields[BED_ITEM_RGB_INDEX];
        parsedData.additionalFields[ITEM_RGB_QUALIFIER_NAME] = itemRgbStr;

        if (false == validateAnnotationColor(itemRgbStr, parsedData.annotColor)) {
            status.incorrectItemRgb = true;
        }
    }

    // Blocks (i.e. exons) parameters
    if (numOfFields > BED_BLOCK_COUNT_INDEX) {
        CHECK_FIELD(BED_BLOCK_COUNT_INDEX);
        QString blockCountStr = fields[BED_BLOCK_COUNT_INDEX];
        QString blockSizesStr;
        QString blockStartsStr;

        // If they are present, then all three value must be present in a line
        if (numOfFields <= BED_BLOCK_SIZES_INDEX) {
            status.incorrectBlocks = true;
        }
        else {
            CHECK_FIELD(BED_BLOCK_SIZES_INDEX);
            blockSizesStr = fields[BED_BLOCK_SIZES_INDEX];

            if (numOfFields <= BED_BLOCK_STARTS_INDEX) {
                status.incorrectBlocks = true;
            }
            else {
                CHECK_FIELD(BED_BLOCK_STARTS_INDEX);
                blockStartsStr = fields[BED_BLOCK_STARTS_INDEX];
            }
        }

        parsedData.additionalFields[BLOCK_COUNT_QUALIFIER_NAME] = blockCountStr;
        parsedData.additionalFields[BLOCK_SIZES_QULAIFIER_NAME] = blockSizesStr;
        parsedData.additionalFields[BLOCK_STARTS_QUALIFIER_NAME] = blockStartsStr;

        // Validate the values
        if (false == validateBlocks(blockCountStr, blockSizesStr, blockStartsStr, parsedData.region)) {
            status.incorrectBlocks = true;
        }
    }

    return parsedData;
}


FormatCheckResult BedFormat::checkRawData(const QByteArray& rawData, const GUrl& /* = GUrl */) const
{
    const char* data = rawData.constData();
    int size = rawData.size();

    bool hasBinaryData = TextUtils::contains(TextUtils::BINARY, data, size);
    if (hasBinaryData) {
        return FormatDetection_NotMatched;
    }

    QString dataStr(rawData);
    QStringList fileLines = dataStr.split("\n");
    BEDLineValidateFlags validationStatus;

    int numToIterate;
    int HUGE_DATA = 65536;
    if (size < HUGE_DATA) {
        numToIterate = fileLines.size(); 
    }
    else {
        // Skip the last line as it can be incomplete
        numToIterate = fileLines.size() - 1;
    }

    bool trackLineDetected = false; // A line that starts with "track" keyword should be present
    int numberOfFieldsPerLine = 0;
    bool firstAnnotLine = true;
    for (int i = 0; i < numToIterate; ++i) {
        if (!fileLines[i].trimmed().isEmpty()) { // e.g. the last line in file can be empty

            QString line = fileLines[i];
            // Skip the header
            if (line.startsWith("browser")) {
                continue;
            }
            if (line.startsWith("track")) {
                trackLineDetected = true;
                continue;
            }

            // Validate other lines
            if (trackLineDetected) {

                // The number of fields per line is detected from the first line
                if (firstAnnotLine) {
                    firstAnnotLine = false;
                    numberOfFieldsPerLine = line.split("\t").count();
                    // There must be at least "chrom", "chromStart" and "chromEnd" fields
                    if (numberOfFieldsPerLine < 3) {
                        return FormatDetection_NotMatched;
                    }
                }

                // Parse a line
                parseAndValidateLine(line, numberOfFieldsPerLine, validationStatus);
            }
        }
    }
    validationStatus.hasTrackLine = trackLineDetected;

    return validationStatus.getFormatDetectionScore( );
}

#define READ_BUFF_SIZE 4096
int readBedLine(QString &buffer, IOAdapter* io, QScopedArrayPointer<char> &charbuff) {
    int len;
    buffer.clear();
    do {
        len = io->readLine(charbuff.data( ), READ_BUFF_SIZE -1);
        charbuff.data( )[len] = '\0';
        buffer.append( QString(charbuff.data( ) ) );
    } while (len == READ_BUFF_SIZE - 1);
    return buffer.length();
}


QList<SharedAnnotationData> BedFormat::getAnnotData( IOAdapter *io, U2OpStatus &os) {
    BedFormat bedFormat( NULL );
    QString annotName = "misc_feature";
    QList<SharedAnnotationData> res;
    const QHash<QString, QList<SharedAnnotationData> > &resHash
        = bedFormat.parseDocument(io, annotName, os);
    foreach ( const QString &seqName, resHash.keys( ) ){
        res.append( resHash.value(seqName ) );
    }
    return res;
}

/**
 * Gets an attribute value in one of the following formats:
 *   1) name="value"
 *   2) name=value[SPACE CHARACTER]
 * If the value is found or the attribute is absent, returns true.
 * If the format of the value is incorrect, returns false.
 */
bool getAttributeValue(const QString& line, const QString& attrName, QString& attrValue)
{
    QString attrStr= attrName + "=";
    int attrIndex = line.indexOf(attrStr);

    if (-1 == attrIndex) {
        // The attribute is not found
        return true;
    }

    int attrBeginIndex = attrIndex + attrStr.size();

    if (line.size() == attrBeginIndex) {
        // Format is incorrect: attribute has name and assignment sign, but not the value
        return false;
    }

    bool parenthesisAreUsed;
    int attrEndIndex;
    if (line[attrBeginIndex] == '\"') {
        parenthesisAreUsed = true;
        attrEndIndex = line.indexOf("\"", attrBeginIndex + 1);
    }
    else {
        parenthesisAreUsed = false;
        attrEndIndex = line.indexOf(QRegExp("\\s"), attrIndex);
        if (-1 == attrEndIndex) {
            attrEndIndex = line.size();
        }
    }

    if (-1 != attrEndIndex) {
        if (parenthesisAreUsed) {
            attrBeginIndex++;
        }
        attrValue = line.mid(attrBeginIndex, attrEndIndex - attrBeginIndex);
        return true;
    }
    else {
        // Format is incorrect: there is no matching end character for the attribute
        return false;
    }
}


/** Get name and description from the track line */
bool parseTrackLine(const QString& trackLine, QString& trackName, QString& trackDescr)
{
    SAFE_POINT(trackLine.startsWith("track "), "Internal error while parsing track header line of a BED file:"
        " the line doesn't starts with 'track'!", false);

    bool attrFormatStatus = true; // Correct, by default

    attrFormatStatus = getAttributeValue(trackLine, "name", trackName);
    if (!attrFormatStatus) {
        return false;
    }

    attrFormatStatus = getAttributeValue(trackLine, "description", trackDescr);

    return attrFormatStatus;
}



QHash<QString, QList<SharedAnnotationData> > BedFormat::parseDocument(
    IOAdapter *io, const QString &defaultAnnotName, U2OpStatus &os )
{
    QString seqName;
    QHash<QString, QList<SharedAnnotationData> > resultHash;
    QList<SharedAnnotationData> result;

    int length;
    QScopedArrayPointer<char> buff( new char[READ_BUFF_SIZE] );
    QString qstrbuf;

    bool fileIsValid = true;
    bool noHeader = false;

    // Parse and validate the header: ignore lines with "browser"
    // Search the 'track' line and get parameters from it
    QString trackName;
    QString trackDescr;

    bool headerLine = true;
    while (headerLine && (length = readBedLine(qstrbuf, io, buff)) > 0) {
        if (qstrbuf.startsWith("#")){ //skip comments
            continue;
        }
        if (qstrbuf.startsWith("browser")) {
            continue;
        }
        else if (qstrbuf.startsWith("track")) {
            if (false == parseTrackLine(qstrbuf, trackName, trackDescr)) {
                fileIsValid = false;
                ioLog.trace(tr("BED parsing error: incorrect format of the 'track' header line!"));
            }
            break; // Stop parsing the header when 'track' line has been detected
        }
        else {
            noHeader = true;
            break;
        }
    }

    // Read other lines
    int lineNumber = 1;
    int numOfFieldsPerLine = 0;
    
    //we have already red the line if there is no header
    if (!noHeader){
        length = readBedLine(qstrbuf, io, buff);
    }

    while (length > 0) {
        
        // Parse and validate the line
        BEDLineValidateFlags validationStatus;

        if (qstrbuf.startsWith("#")){//skip comments
            os.setProgress(io->getProgress());
            length = readBedLine(qstrbuf, io, buff);
            continue;
        }
        if (1 == lineNumber) {
            numOfFieldsPerLine = qstrbuf.split("\t").count();
            // "3" as there must be at least "chrom", "chromStart" and "chromEnd" fields
            if (numOfFieldsPerLine < 3) {
                os.setError(tr("BED parsing error: unexpected number of fields in the first annotations line!"));
            }
        }
        BedLineData bedLineData = parseAndValidateLine(qstrbuf, numOfFieldsPerLine, validationStatus);

        // Check that an annotation can be created
        if (true == validationStatus.incorrectNumberOfFields) {
            os.setError(tr("BED parsing error: incorrect number of fields at line %1!").arg(lineNumber));
            addToResults(resultHash, result, seqName);
            return resultHash;
        }

        if (true == validationStatus.emptyFields) {
            os.setError(tr("BED parsing error: a field at line %1 is empty!").arg(lineNumber));
            addToResults(resultHash, result, seqName);
            return resultHash;
        }

        if (true == validationStatus.incorrectCoordinates) {
            os.setError(tr("BED parsing error: incorrect coordinates at line %1!").arg(lineNumber));
            addToResults(resultHash, result, seqName);
            return resultHash;
        }

        // If file is invalid, but can be parsed an error is written to the log,
        // all details are written to the trace log.
        if (validationStatus.isFileInvalid() && !noHeader) {
            fileIsValid = false;
        }

        if (!seqName.isEmpty()) {
            if (bedLineData.seqName != seqName) {
                addToResults(resultHash, result, seqName);
                seqName = bedLineData.seqName;
            }
        }
        else {
           seqName = bedLineData.seqName;
       }

        // Create the annotation
        SharedAnnotationData annotData(new AnnotationData());
        annotData->name = defaultAnnotName;
        annotData->location->regions << bedLineData.region;

        // Add qualifiers
        foreach (QString qualifierName, bedLineData.additionalFields.keys()) {
            annotData->qualifiers.push_back(
                U2Qualifier(qualifierName, bedLineData.additionalFields.value(qualifierName)));
        }

        // Add a qualifier with the sequence name
        annotData->qualifiers.push_back(U2Qualifier(CHROM_QUALIFIER_NAME, bedLineData.seqName));

        // Additionally, verify the strand information
        if (bedLineData.additionalFields.keys().contains(STRAND_QUALIFIER_NAME)) {
            if ("-" == bedLineData.additionalFields.value(STRAND_QUALIFIER_NAME)) {
                annotData->setStrand(U2Strand::Complementary);
            }
        }

        // Add track information to the qualifiers
        if (!trackName.isEmpty()) {
            annotData->qualifiers.push_back(U2Qualifier(TRACK_NAME_QUALIFIER_NAME, trackName));
        }

        if (!trackDescr.isEmpty()) {
            annotData->qualifiers.push_back(U2Qualifier(TRACK_DESCR_QUALIFIER_NAME, trackDescr));
        }

        // If there were some errors during parsing the output, write it to the log
        if (true == validationStatus.incorrectScore) {
            ioLog.trace(tr("BED parsing error: incorrect score value '%1'"
                " at line %2!").arg(bedLineData.additionalFields[SCORE_QUALIFIER_NAME]).arg(lineNumber));
        }

        if (true == validationStatus.incorrectStrand) {
            ioLog.trace(tr("BED parsing error: incorrect strand value '%1'"
                " at line %2!").arg(bedLineData.additionalFields[STRAND_QUALIFIER_NAME]).arg(lineNumber));
        }

        if (true == validationStatus.incorrectThickCoordinates) {
            ioLog.trace(tr("BED parsing error: incorrect thick coordinates at line %1!").arg(lineNumber));
        }

        if (true == validationStatus.incorrectItemRgb) {
            ioLog.trace(tr("BED parsing error: incorrect itemRgb value '%1'"
                " at line %2!").arg(bedLineData.additionalFields[ITEM_RGB_QUALIFIER_NAME]).arg(lineNumber));
        }

        if (true == validationStatus.incorrectBlocks) {
            ioLog.trace(tr("BED parsing error: incorrect value of the block parameters"
                " at line %1!").arg(lineNumber));
        }

        // Append the result
        result.append(annotData);

        // Move to the next line
        lineNumber++;

        os.setProgress(io->getProgress());

        length = readBedLine(qstrbuf, io, buff);
    }

    if (false == fileIsValid) {
        ioLog.info("BED parsing warning: one or more errors occurred while parsing the input file,"
            " see TRACE log for details!");
    }

    addToResults(resultHash, result, seqName);
    return resultHash;
}

void BedFormat::storeDocument(Document* doc, IOAdapter* io, U2OpStatus& os)
{
    SAFE_POINT(NULL != doc, "Internal error: NULL document was provided to BEDFormat::storeDocument!",);
    SAFE_POINT(NULL != io, "Internal error: NULL IO adapter was provided to BEDFormat::storeDocument!",);

    ioLog.trace(tr("Starting BED saving: '%1'").arg(doc->getURLString()));

    bool noErrorsDuringStoring = true;
    QList<GObject*> annotTables = doc->findGObjectByType(GObjectTypes::ANNOTATION_TABLE);

    QByteArray lineData;

    int fieldsNumberPerLine = 0;
    bool firstLine = true;

    foreach (GObject* annotTable, annotTables) {
        QList<Annotation> annotationsList =
            ( qobject_cast<AnnotationTableObject *>( annotTable ) )->getAnnotations( );

        foreach (const Annotation &annot, annotationsList) {
            QString annotName = annot.getName( );
            if (annotName == U1AnnotationUtils::lowerCaseAnnotationName ||
                annotName == U1AnnotationUtils::upperCaseAnnotationName)
            {
                continue;
            }

            QStringList lineFields;
            QVector<U2Region> annotRegions = annot.getRegions( );
            QVector<U2Qualifier> annotQualifiers = annot.getQualifiers( );

            QString chromName = annot.findFirstQualifierValue( CHROM_QUALIFIER_NAME );
            if (chromName.isEmpty()) {
                ioLog.trace(tr("BED saving error: can't save an annotation to a BED file"
                    " - the annotation doesn't have the 'chrom' qualifier!"));
                noErrorsDuringStoring = false;
                continue;
            }
            else {
                lineFields << chromName;
            }

            foreach (const U2Region& region, annotRegions) {
                // chromStart and chromEnd
                lineFields << QString::number(region.startPos);
                lineFields << QString::number(region.endPos());

                QString nameQualValue = annot.findFirstQualifierValue(ANNOT_QUALIFIER_NAME);
                QString scoreQualValue = annot.findFirstQualifierValue(SCORE_QUALIFIER_NAME);
                QString strandQualValue = annot.findFirstQualifierValue(STRAND_QUALIFIER_NAME);
                QString thickStartQualValue = annot.findFirstQualifierValue(THICK_START_QUALIFIER_NAME);
                QString thickEndQualValue = annot.findFirstQualifierValue(THICK_END_QUALIFIER_NAME);
                QString itemRgbQualValue = annot.findFirstQualifierValue(ITEM_RGB_QUALIFIER_NAME);
                QString blockCountQualValue = annot.findFirstQualifierValue(BLOCK_COUNT_QUALIFIER_NAME);
                QString blockSizesQualValue = annot.findFirstQualifierValue(BLOCK_SIZES_QULAIFIER_NAME);
                QString blockStartsQualValue = annot.findFirstQualifierValue(BLOCK_STARTS_QUALIFIER_NAME);

                // Write the header, detect the number of optional fields from the first annotation
                if (firstLine) {
                    firstLine = false;
                    // Header
                    QString trackNameQualValue = annot.findFirstQualifierValue(TRACK_NAME_QUALIFIER_NAME);
                    QString trackDescrQualValue = annot.findFirstQualifierValue(TRACK_DESCR_QUALIFIER_NAME);

                    if ( ! ( trackNameQualValue.isEmpty( ) || trackDescrQualValue.isEmpty( ) ) ) {
                        QString headerStr = QString( "track name=\"%1\" description=\"%2\"\n" )
                            .arg( trackNameQualValue ).arg( trackDescrQualValue );
                        QByteArray header = headerStr.toLatin1();
                        qint64 len = io->writeBlock(header);
                        if (len != header.size()) {
                            os.setError(L10N::errorWritingFile(doc->getURLString()));
                            return;
                        }
                    }

                    // Number of optional fields
                    // Note that the order of the optional fields is binding:
                    // lower-numbered fields must always be populated
                    // if higher-numbered fields are used
                    if (nameQualValue.isEmpty()) {
                        fieldsNumberPerLine = 3; // No default value, skip all optional fields
                    }
                    else {
                        fieldsNumberPerLine = 4;

                        // If score and strand qualifiers are not present, but further qualifiers are present,
                        // the line is filled with values "0" and "+" (or "-" depending on the annotation strand) at these positions
                        if (!scoreQualValue.isEmpty()) {
                            fieldsNumberPerLine = 5;
                        }
                        if (!strandQualValue.isEmpty()) {
                            fieldsNumberPerLine = 6;
                        }

                        // If thick coordinates qualifier are not present, the annotation
                        // region's coordinates are used
                        if (!thickStartQualValue.isEmpty()) {
                            if (!thickEndQualValue.isEmpty()) {
                                fieldsNumberPerLine = 8;
                            }
                            else {
                                ioLog.trace(tr("BED saving error: incorrect thick coordinates"
                                    " in the first annotation!"));
                                noErrorsDuringStoring = false;
                            }
                        }

                        // Red color (255, 0, 0) is used by default for itemRgb
                        if (!itemRgbQualValue.isEmpty()) {
                            fieldsNumberPerLine = 9;
                        }

                        if (!blockCountQualValue.isEmpty()) {
                            if (blockStartsQualValue.isEmpty() ||
                                blockSizesQualValue.isEmpty()) {
                                    ioLog.trace(tr("BED saving error: incorrect block fields"
                                        " in the first annotation!"));
                                    noErrorsDuringStoring = false;
                            }
                            else {
                                fieldsNumberPerLine = 12;
                            }
                        }
                    }
                    ioLog.trace(tr("BED saving: detected %1 fields per line"
                        " for file '%2'").arg(fieldsNumberPerLine).arg(doc->getURLString()));
                }

                // Append the required number of fields to the line
                if (fieldsNumberPerLine >= 4) {
                    if (nameQualValue.isEmpty()) {
                        ioLog.trace(tr("BED saving error: an annotation is expected to have '%1'"
                            " qualifier, but it is absent! Skipping the annotation.").arg(ANNOT_QUALIFIER_NAME));
                        noErrorsDuringStoring = false;
                        continue;
                    }
                    else {
                        lineFields << nameQualValue;
                    }
                }

                if (fieldsNumberPerLine >= 5) {
                    if (scoreQualValue.isEmpty()) {
                        lineFields << "0";
                    }
                    else {
                        lineFields << scoreQualValue;
                    }
                }

                if (fieldsNumberPerLine >= 6) {
                    if (strandQualValue.isEmpty()) {
                        U2Strand strand = annot.getStrand();
                        if (strand == U2Strand::Complementary) {
                            lineFields << "-";
                        }
                        else {
                            lineFields << "+";
                        }
                    }
                    else {
                        lineFields << strandQualValue;
                    }
                }

                if (fieldsNumberPerLine >= 8) {
                    if (thickStartQualValue.isEmpty() ||
                        thickEndQualValue.isEmpty()) {
                            // Write chromStart and chromEnd coordinates
                            lineFields << QString::number(region.startPos);
                            lineFields << QString::number(region.endPos());
                    }
                    else {
                        lineFields << thickStartQualValue;
                        lineFields << thickEndQualValue;
                    }
                }

                if (fieldsNumberPerLine >= 9) {
                    if (itemRgbQualValue.isEmpty()) {
                        lineFields << "255, 0, 0";
                    }
                    else {
                        lineFields << itemRgbQualValue;
                    }
                }

                if (fieldsNumberPerLine >= 12) {
                    if (blockCountQualValue.isEmpty() ||
                        blockStartsQualValue.isEmpty() ||
                        blockSizesQualValue.isEmpty())
                    {
                        ioLog.trace(tr("BED saving error: an annotation is expected to have the block"
                            " qualifiers! Skipping the annotation."));
                        noErrorsDuringStoring = false;
                        continue;
                    }
                    else {
                        lineFields << blockCountQualValue;
                        lineFields << blockSizesQualValue;
                        lineFields << blockStartsQualValue;
                    }
                }

                // Write the line
                lineData = lineFields.join("\t").toLatin1() + "\n";
                qint64 len = io->writeBlock(lineData);
                if (len != lineData.size()) {
                    os.setError(L10N::errorWritingFile(doc->getURLString()));
                    return;
                }
            }
        }
    }

    if (!noErrorsDuringStoring) {
        ioLog.error(tr("BED saving error: one or more errors occurred while saving file '%1',"
            " see TRACE log for details!").arg(doc->getURLString()));
    }

    ioLog.trace(tr("Finished BED saving: '%1'").arg(doc->getURLString()));
}

void BedFormat::addToResults( QHash<QString, QList<SharedAnnotationData> > & resHash, QList<SharedAnnotationData>& result, const QString& seqName ){
    QHash<QString, QList<SharedAnnotationData> >::iterator i = resHash.find(seqName);
    if (i != resHash.end()) {
        i.value().append(result);
    }else{
        resHash.insert(seqName, result);
    }
    result.clear();
}

} // namespace U2
