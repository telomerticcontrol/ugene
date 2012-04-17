/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2012 UniPro <ugene@unipro.ru>
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

#include "CollocationWorker.h"

#include <U2Lang/IntegralBusModel.h>
#include <U2Lang/WorkflowEnv.h>
#include <U2Lang/ActorPrototypeRegistry.h>
#include <U2Core/QVariantUtils.h>
#include <U2Lang/BaseTypes.h>
#include <U2Lang/BaseSlots.h>
#include <U2Lang/BasePorts.h>
#include <U2Lang/BaseActorCategories.h>
#include <U2Designer/DelegateEditors.h>
#include <U2Lang/CoreLibConstants.h>

#include <U2Core/AnnotationData.h>
#include <U2Core/DNASequence.h>
#include <U2Core/DNASequenceObject.h>
#include <U2Core/DNATranslation.h>
#include <U2Core/DNAAlphabet.h>
#include <U2Core/AppContext.h>
#include <U2Core/Log.h>

#include "CollocationsSearchAlgorithm.h"
#include "CollocationsDialogController.h"

/* TRANSLATOR U2::LocalWorkflow::CollocationWorker */

namespace U2 {
namespace LocalWorkflow {

static const QString NAME_ATTR("result-name");
static const QString ANN_ATTR("annotations");
static const QString LEN_ATTR("region-size");
static const QString FIT_ATTR("must-fit");

const QString SEQ_SLOT = BaseSlots::DNA_SEQUENCE_SLOT().getId();
const QString FEATURE_TABLE_SLOT = BaseSlots::ANNOTATION_TABLE_SLOT().getId();

const QString CollocationWorkerFactory::ACTOR_ID("collocated-annotation-search");

class CollocationValidator : public ConfigurationValidator {
public:
    virtual bool validate(const Configuration* cfg, QStringList& output) const {
        QString annotations = cfg->getParameter(ANN_ATTR)->getAttributeValueWithoutScript<QString>();
        QSet<QString> names = QSet<QString>::fromList(annotations.split(QRegExp("\\W+"), QString::SkipEmptyParts));
        if (names.size() < 2) {
            output.append(CollocationWorker::tr("At least 2 annotations are required for collocation search."));
            return false;
        }
        return true;
    }
};

void CollocationWorkerFactory::init() {

    QMap<Descriptor, DataTypePtr> m;
    {
        //Descriptor sd(DNA_SEQUENCE_SLOT, CollocationWorker::tr("sequence-slot"), CollocationWorker::tr("seq-slot-doc"));
        //Descriptor fd(ANNOTATION_TABLE_SLOT, CollocationWorker::tr("feature-table"), CollocationWorker::tr("feature-table-doc"));
        m[BaseSlots::DNA_SEQUENCE_SLOT()] = BaseTypes::DNA_SEQUENCE_TYPE();
        m[BaseSlots::ANNOTATION_TABLE_SLOT()] = BaseTypes::ANNOTATION_TABLE_LIST_TYPE();
    }
    DataTypePtr inSet(new MapDataType(Descriptor("regioned.sequence"), m));
    DataTypeRegistry* dr = WorkflowEnv::getDataTypeRegistry();
    assert(dr);
    dr->registerEntry(inSet);

    QList<PortDescriptor*> p; QList<Attribute*> a;
    p << new PortDescriptor(Descriptor(BasePorts::IN_SEQ_PORT_ID(), CollocationWorker::tr("Input data"), 
        CollocationWorker::tr("An input sequence and a set of annotations to search in.")), inSet, true /*input*/);
    QMap<Descriptor, DataTypePtr> outM; 
    outM[BaseSlots::ANNOTATION_TABLE_SLOT()] = BaseTypes::ANNOTATION_TABLE_TYPE();
    p << new PortDescriptor(Descriptor(BasePorts::OUT_ANNOTATIONS_PORT_ID(), 
        CollocationWorker::tr("Group annotations"), CollocationWorker::tr("Annotated regions containing found collocations.")),
        DataTypePtr(new MapDataType(Descriptor("collocation.annotations"), outM)), false /*input*/, true/*multi*/);
    
    {
        //U2Region             searchRegion;
        Descriptor nd(NAME_ATTR, CollocationWorker::tr("Result annotation"), CollocationWorker::tr("Name of the result annotations to mark found collocations"));
        Descriptor ad(ANN_ATTR, CollocationWorker::tr("Group of annotations"), CollocationWorker::tr("A list of annotation names to search. Found regions will contain all the named annotations."));
        Descriptor ld(LEN_ATTR, CollocationWorker::tr("Region size"), CollocationWorker::tr("Effectively this is the maximum allowed distance between the interesting annotations in a group"));
        Descriptor fd(FIT_ATTR, CollocationWorker::tr("Must fit into region"), CollocationWorker::tr("Whether the interesting annotations should entirely fit into the specified region to form a group"));
        a << new Attribute(nd, BaseTypes::STRING_TYPE(), true, QVariant("misc_feature"));
        a << new Attribute(ad, BaseTypes::STRING_TYPE(), true);
        a << new Attribute(ld, BaseTypes::NUM_TYPE(), false, QVariant(1000));
        a << new Attribute(fd, BaseTypes::BOOL_TYPE(), false, QVariant(false));
    }

    Descriptor desc(ACTOR_ID, CollocationWorker::tr("Collocation Search"), 
        CollocationWorker::tr("Finds groups of specified annotations in each supplied set of annotations, stores found regions as annotations."));
    ActorPrototype* proto = new IntegralBusActorPrototype(desc, p, a);
    proto->addSlotRelation(BasePorts::IN_SEQ_PORT_ID(), BaseSlots::DNA_SEQUENCE_SLOT().getId(),
        BasePorts::OUT_ANNOTATIONS_PORT_ID(), BaseSlots::ANNOTATION_TABLE_SLOT().getId());
    QMap<QString, PropertyDelegate*> delegates;    
    
    QVariantMap lenMap; lenMap["minimum"] = QVariant(0); lenMap["maximum"] = QVariant(INT_MAX);
    delegates[LEN_ATTR] = new SpinBoxDelegate(lenMap);
       
    proto->setEditor(new DelegateEditor(delegates));
    proto->setValidator(new CollocationValidator());
    proto->setIconPath(":annotator/images/regions.png");
    proto->setPrompter(new CollocationPrompter());
    WorkflowEnv::getProtoRegistry()->registerProto(BaseActorCategories::CATEGORY_BASIC(), proto);

    DomainFactory* localDomain = WorkflowEnv::getDomainRegistry()->getById(LocalDomainFactory::ID);
    localDomain->registerEntry(new CollocationWorkerFactory());
}

QString CollocationPrompter::composeRichDoc() {
    IntegralBusPort* input = qobject_cast<IntegralBusPort*>(target->getPort(BasePorts::IN_SEQ_PORT_ID()));
    Actor* seqProducer = input->getProducer(SEQ_SLOT);
    QString seqName = seqProducer ? tr(" sequence from <u>%1</u>").arg(seqProducer->getLabel()) : "";
    QString annName = getProducers(BasePorts::IN_SEQ_PORT_ID(), FEATURE_TABLE_SLOT);
    if (!annName.isEmpty()) {
        annName = tr(" set of annotations from <u>%1</u>").arg(annName);
    }

    QString data;
    if (seqName.isEmpty() && annName.isEmpty()) {
        //return "<font color='red'>"+tr("unset")+"</font>";
    } else if (!seqName.isEmpty() && !annName.isEmpty()) {
        data = tr("For each %1 and %2,").arg(seqName).arg(annName);
    } else {
        data = tr("For each %1%2,").arg(seqName).arg(annName);
    }

    QString annotations;
    QStringList names = annotations.split(QRegExp("\\W+"), QString::SkipEmptyParts).toSet().toList();
    annotations = names.join(", ");
    if (annotations.isEmpty()) {
        annotations = getRequiredParam(ANN_ATTR);
    }
    annotations = getHyperlink(ANN_ATTR, annotations);

    int distance = getParameter(LEN_ATTR).toInt();
    bool mode = getParameter(FIT_ATTR).toBool();

    QString extra;
    if (mode) {
        extra = tr(" Annotations themselves may not span beyond the region.");
    }

    QString resultName = getHyperlink(NAME_ATTR, getRequiredParam(NAME_ATTR));
    QString doc = tr("%1 look if <u>%2</u> annotations appear collocated within same region of length <u>%3</u>.%4"
        "<br>Output the list of found regions annotated as <u>%5</u>.")
        .arg(data) //sequence from Read Fasta 1
        .arg(annotations)
        .arg(getHyperlink(LEN_ATTR, distance))
        .arg(extra)
        .arg(resultName);

    return doc;
}

void CollocationWorker::init() {
    input = ports.value(BasePorts::IN_SEQ_PORT_ID());
    output = ports.value(BasePorts::OUT_ANNOTATIONS_PORT_ID());
}

Task* CollocationWorker::tick() {
    if (input->hasMessage()) {
        Message inputMessage = getMessageAndSetupScriptValues(input);
        if (inputMessage.isEmpty()) {
            output->transit();
            return NULL;
        }
        cfg.distance = actor->getParameter(LEN_ATTR)->getAttributeValue<int>(context);
        cfg.st = actor->getParameter(FIT_ATTR)->getAttributeValue<bool>(context) ? 
            CollocationsAlgorithm::NormalSearch : CollocationsAlgorithm::PartialSearch;
        resultName = actor->getParameter(NAME_ATTR)->getAttributeValue<QString>(context);
        QString annotations = actor->getParameter(ANN_ATTR)->getAttributeValue<QString>(context);
        names = QSet<QString>::fromList(annotations.split(QRegExp("\\W+"), QString::SkipEmptyParts));
        QVariantMap qm = inputMessage.getData().toMap();

        U2DataId seqId = qm.value(BaseSlots::DNA_SEQUENCE_SLOT().getId()).value<U2DataId>();
        std::auto_ptr<U2SequenceObject> seqObj(StorageUtils::getSequenceObject(context->getDataStorage(), seqId));
        if (NULL == seqObj.get()) {
            return NULL;
        }
        DNASequence seq = seqObj->getWholeSequence();
        
        QList<SharedAnnotationData> atl = QVariantUtils::var2ftl(qm.value(FEATURE_TABLE_SLOT).toList());
        if (!seq.isNull() && !atl.isEmpty()) {
            cfg.searchRegion.length = seq.length();
            Task* t = new CollocationSearchTask(atl, names, cfg);
            connect(t, SIGNAL(si_stateChanged()), SLOT(sl_taskFinished()));
            return t;
        } else {
            // void tick
            output->put(Message(BaseTypes::ANNOTATION_TABLE_TYPE(), QVariant()));
            if (input->isEnded()) {
                output->setEnded();
            }
            return NULL;
        }
    } else if (input->isEnded()) {
        setDone();
        output->setEnded();
    }
    return NULL;
}

void CollocationWorker::sl_taskFinished() {
    CollocationSearchTask* t = qobject_cast<CollocationSearchTask*>(sender());
    if (t->getState() != Task::State_Finished) return;
    QVector<U2Region> res = t->popResults();
    if (output) {
        QList<SharedAnnotationData> list;
        foreach(U2Region r, res) {
            SharedAnnotationData data; data = new AnnotationData();
            data->location->regions.append(r);
            data->setStrand(U2Strand::Direct);
            data->name = resultName;
            list.append(data);
        }

        QVariant v = qVariantFromValue<QList<SharedAnnotationData> >(list);
        output->put(Message(BaseTypes::ANNOTATION_TABLE_TYPE(), v));
        algoLog.info(tr("Found %1 collocations").arg(res.size()));
    }
}

} //namespace LocalWorkflow
} //namespace U2
