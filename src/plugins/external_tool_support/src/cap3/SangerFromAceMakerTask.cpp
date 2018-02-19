/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2018 UniPro <ugene@unipro.ru>
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

#include "SangerFromAceMakerTask.h"

#include <U2Core/AppContext.h>
#include <U2Core/BaseDocumentFormats.h>
#include <U2Core/ChromatogramUtils.h>
#include <U2Core/DNASequenceObject.h>
#include <U2Core/DocumentModel.h>
#include <U2Core/GObjectReference.h>
#include <U2Core/GObjectRelationRoles.h>
#include <U2Core/IOAdapterUtils.h>
#include <U2Core/MsaDbiUtils.h>
#include <U2Core/MultipleChromatogramAlignmentImporter.h>
#include <U2Core/MultipleChromatogramAlignmentObject.h>
#include <U2Core/SaveDocumentTask.h>
#include <U2Core/U2AlphabetUtils.h>
#include <U2Core/U2AttributeDbi.h>
#include <U2Core/U2AttributeUtils.h>
#include <U2Core/U2ObjectDbi.h>
#include <U2Core/U2OpStatusUtils.h>
#include <U2Core/U2SafePoints.h>
#include <U2Core/U2SequenceUtils.h>

#include <U2Formats/AceImportUtils.h>


namespace U2 {

SangerFromAceMakerTask::SangerFromAceMakerTask(const GUrl& aceFile,
                                               const QMap<QString, DNAChromatogram>& chromMap,
                                               const GUrl& outputFile)
    : Task("Make ugenedb from ace and a bunch of abi files", TaskFlag_None),
      chromMap(chromMap),
      aceFile(aceFile),
      outputFile(outputFile) {

}

void SangerFromAceMakerTask::run() {
    QScopedPointer<IOAdapter> ioAdapter;
    IOAdapterFactory *factory = AppContext::getIOAdapterRegistry()->getIOAdapterFactoryById(IOAdapterUtils::url2io(aceFile));
    SAFE_POINT_EXT(factory, setError(tr("IOAdapterFactory is NULL")), );
    ioAdapter.reset(factory->createIOAdapter());

    if (!ioAdapter->open(aceFile, IOAdapterMode_Read)) {
        setError(tr("Can't open file '%1'").arg(aceFile.getURLString()));
        return;
    }

    QScopedPointer<AceReader> aceReader;
    U2OpStatusChildImpl os(&stateInfo, U2OpStatusMapping(0, 50));
    aceReader.reset(new AceReader(*ioAdapter, os));
    CHECK_OP(os, );

    QScopedPointer<AceIterator> iterator;
    iterator.reset(new AceIterator(*aceReader, stateInfo));

    // TODO: make a few objects in case of several contigs
    MultipleChromatogramAlignment alignment;
    DNASequence refSequence;
    while (iterator->hasNext()) {
        Assembly aceAssembly = iterator->next();
        CHECK_OP(stateInfo, );

        Assembly::Sequence aceReference = aceAssembly.getReference();
        alignment->setName(aceAssembly.getName());
        refSequence = DNASequence("contig", aceReference.data);

        QList<Assembly::Sequence> reads = aceAssembly.getInitialSequences();
        int index = 0;
        foreach (Assembly::Sequence read, reads) {
            QString seqName(read.name);
            if (chromMap.contains(seqName)) {
                QByteArray seq;
                U2MsaRowGapModel gaps;
                // HACK: see lines 387-389 in AceFormat.cpp
                read.data.replace('*',U2Msa::GAP_CHAR);
                read.data.replace('N',U2Msa::GAP_CHAR);
                read.data.replace('X',U2Msa::GAP_CHAR);
                MaDbiUtils::splitBytesToCharsAndGaps(read.data, seq, gaps);
                MsaRowUtils::addOffsetToGapModel(gaps, read.offset);

                if (read.isComplemented) {
                    DNAChromatogram chromatogram = ChromatogramUtils::reverseComplement(chromMap[seqName]);
                    alignment->addRow(seqName, chromatogram, DNASequence(seqName, seq), gaps, stateInfo, read.isComplemented);
                } else {
                    alignment->addRow(seqName, chromMap[seqName], DNASequence(seqName, seq), gaps, stateInfo, read.isComplemented);
                }
                CHECK_OP(stateInfo, );
                index++;
            }
        }
    }
    // TODO: detect the alphabet by contig or/and reads
    alignment->setAlphabet(U2AlphabetUtils::getById(BaseDNAAlphabetIds::NUCL_DNA_EXTENDED()));

    U2DbiRef ref(DEFAULT_DBI_ID, outputFile.getURLString());
    DbiConnection con(ref, true, stateInfo);
    Q_UNUSED(con);
    CHECK_OP(stateInfo, );

    MultipleChromatogramAlignmentObject* mcaObject = MultipleChromatogramAlignmentImporter::createAlignment(stateInfo, ref, U2ObjectDbi::ROOT_FOLDER, alignment);
    CHECK_OP(stateInfo, );
    SAFE_POINT_EXT(NULL != mcaObject, setError("Result MCA object is NULL"), );

    U2EntityRef seqEntityRef = U2SequenceUtils::import(stateInfo, ref, U2ObjectDbi::ROOT_FOLDER,
                                                       refSequence, BaseDNAAlphabetIds::NUCL_DNA_EXTENDED());
    CHECK_OP(stateInfo, );

    U2SequenceObject *referenceSequenceObject = new U2SequenceObject("contig", seqEntityRef);
    SAFE_POINT_EXT(NULL != referenceSequenceObject, setError("Result reference sequence object is NULL"), );

    U2ByteArrayAttribute attribute;
    U2Object obj;
    obj.dbiId = ref.dbiId;
    obj.id = mcaObject->getEntityRef().entityId;
    obj.version = mcaObject->getModificationVersion();
    U2AttributeUtils::init(attribute, obj, MultipleChromatogramAlignmentObject::MCAOBJECT_REFERENCE);
    attribute.value = referenceSequenceObject->getEntityRef().entityId;
    con.dbi->getAttributeDbi()->createByteArrayAttribute(attribute, stateInfo);
    CHECK_OP(stateInfo, );
}

} // namespace
