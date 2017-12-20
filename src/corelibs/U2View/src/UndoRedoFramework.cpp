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

#include <U2Core/AnnotationTableObject.h>
#include <U2Core/DNASequenceObject.h>
#include <U2Core/MultipleSequenceAlignmentObject.h>
#include <U2Core/U2DbiUtils.h>
#include <U2Core/U2ObjectDbi.h>
#include <U2Core/U2OpStatusUtils.h>
#include <U2Core/U2SafePoints.h>

#include <U2Gui/GUIUtils.h>

#include "UndoRedoFramework.h"
#include "ov_msa/MSACollapsibleModel.h"

namespace U2 {


UndoRedoFramework::UndoRedoFramework(QObject *p, GObject *_obj)
    : QObject(p),
      obj(_obj),
      stateComplete(true),
      undoStepsAvailable(0),
      redoStepsAvailable(0)
{
    SAFE_POINT(obj != NULL, "NULL GObject!", );

    undoAction = new QAction(this);
    undoAction->setText(tr("Undo"));
    undoAction->setIcon(QIcon(":core/images/undo.png"));
    undoAction->setShortcut(QKeySequence::Undo);
    GUIUtils::updateActionToolTip(undoAction);

    redoAction = new QAction(this);
    redoAction->setText(tr("Redo"));
    redoAction->setIcon(QIcon(":core/images/redo.png"));
    redoAction->setShortcut(QKeySequence::Redo);
    GUIUtils::updateActionToolTip(redoAction);

    checkUndoRedoEnabled();

    connect(obj, SIGNAL(si_lockedStateChanged()), SLOT(sl_lockedStateChanged()));
    connect(undoAction, SIGNAL(triggered()), this, SLOT(sl_undo()));
    connect(redoAction, SIGNAL(triggered()), this, SLOT(sl_redo()));
}

void UndoRedoFramework::sl_completeStateChanged(bool _stateComplete) {
    stateComplete = _stateComplete;
}

void UndoRedoFramework::sl_lockedStateChanged() {
    checkUndoRedoEnabled();
}

void UndoRedoFramework::checkUndoRedoEnabled() {
    SAFE_POINT(obj != NULL, "NULL MSA Object!", );

    if (obj->isStateLocked() || !stateComplete) {
        undoAction->setEnabled(false);
        redoAction->setEnabled(false);
        return;
    }

    U2OpStatus2Log os;
    DbiConnection con(obj->getEntityRef().dbiRef, os);
    SAFE_POINT_OP(os, );

    U2ObjectDbi* objDbi = con.dbi->getObjectDbi();
    SAFE_POINT(NULL != objDbi, "NULL Object Dbi!", );

    bool enableUndo = objDbi->canUndo(obj->getEntityRef().entityId, os);
    SAFE_POINT_OP(os, );
    bool enableRedo = objDbi->canRedo(obj->getEntityRef().entityId, os);
    SAFE_POINT_OP(os, );

    undoAction->setEnabled(enableUndo);
    redoAction->setEnabled(enableRedo);
}

void UndoRedoFramework::sl_undo() {
    SAFE_POINT(obj != NULL, "NULL MSA Object!", );

    U2OpStatus2Log os;
    U2EntityRef entityRef =  obj->getEntityRef();

    assert(stateComplete);
    assert(!obj->isStateLocked());

    DbiConnection con(entityRef.dbiRef, os);
    SAFE_POINT_OP(os, );

    U2ObjectDbi* objDbi = con.dbi->getObjectDbi();
    SAFE_POINT(NULL != objDbi, "NULL Object Dbi!", );

    QHash<QString, QString> meta = objDbi->undo(entityRef.entityId, os);
    SAFE_POINT_OP(os, );

    updateObject(meta, MaModificationType_Undo);
}

void UndoRedoFramework::sl_redo() {
    SAFE_POINT(obj != NULL, "NULL MSA Object!", );

    U2OpStatus2Log os;
    U2EntityRef entityRef = obj->getEntityRef();

    assert(stateComplete);
    assert(!obj->isStateLocked());

    DbiConnection con(entityRef.dbiRef, os);
    SAFE_POINT_OP(os, );

    U2ObjectDbi* objDbi = con.dbi->getObjectDbi();
    SAFE_POINT(NULL != objDbi, "NULL Object Dbi!", );

    QHash<QString, QString> meta = objDbi->redo(entityRef.entityId, os);
    SAFE_POINT_OP(os, );

    updateObject(meta, MaModificationType_Redo);
}

SequenceUndoRedoFramework::SequenceUndoRedoFramework(QObject *p, U2SequenceObject *seqObj, QList<AnnotationTableObject*> annTableList)
    : UndoRedoFramework(p, seqObj),
      annTableList(annTableList)
{
    connect(seqObj, SIGNAL(si_sequenceChanged()), SLOT(sl_sequenceChanged()));
    foreach(AnnotationTableObject* annTableObj, annTableList) {
        connect(annTableObj, SIGNAL(si_onAnnotationModified(AnnotationModification)), SLOT(sl_lockedStateChanged()));
        connect(annTableObj, SIGNAL(si_modifiedStateChanged()), SLOT(sl_lockedStateChanged()));
        connect(annTableObj, SIGNAL(si_nameChanged(QString)), SLOT(sl_lockedStateChanged()));
        connect(annTableObj, SIGNAL(si_onAnnotationsRemoved(QList<Annotation*>)), SLOT(sl_lockedStateChanged()));
        connect(annTableObj, SIGNAL(si_onAnnotationsAdded(QList<Annotation*>)), SLOT(sl_lockedStateChanged()));
    }

    // TODO_SVEDIT: seq object should notify about changes
    connect(undoAction, SIGNAL(triggered(bool)), SLOT(sl_sequenceChanged()));
    connect(redoAction, SIGNAL(triggered(bool)), SLOT(sl_sequenceChanged()));

    // get adv context
//    v
}

void SequenceUndoRedoFramework::sl_sequenceChanged() {
    checkUndoRedoEnabled();
}

void SequenceUndoRedoFramework::sl_annTableAdded(AnnotationTableObject* annTableObj) {
    connect(annTableObj, SIGNAL(si_onAnnotationModified(AnnotationModification)), SLOT(sl_lockedStateChanged()));
    connect(annTableObj, SIGNAL(si_modifiedStateChanged()), SLOT(sl_lockedStateChanged()));
    connect(annTableObj, SIGNAL(si_nameChanged(QString)), SLOT(sl_lockedStateChanged()));
    connect(annTableObj, SIGNAL(si_onAnnotationsRemoved(QList<Annotation*>)), SLOT(sl_lockedStateChanged()));
    connect(annTableObj, SIGNAL(si_onAnnotationsAdded(QList<Annotation*>)), SLOT(sl_lockedStateChanged()));

    U2OpStatusImpl os;
    annTableObj->setTrackMod(os, TrackOnUpdate);
    SAFE_POINT_OP(os, );

    annTableList.append(annTableObj);
}

void SequenceUndoRedoFramework::sl_annTableRemoved(AnnotationTableObject* annTableObj) {
    disconnect(annTableObj, SIGNAL(si_onAnnotationModified(AnnotationModification)), this, SLOT(sl_lockedStateChanged()));
    disconnect(annTableObj, SIGNAL(si_modifiedStateChanged()), this, SLOT(sl_lockedStateChanged()));
    disconnect(annTableObj, SIGNAL(si_nameChanged(QString)), this, SLOT(sl_lockedStateChanged()));
    disconnect(annTableObj, SIGNAL(si_onAnnotationsRemoved(QList<Annotation*>)), this, SLOT(sl_lockedStateChanged()));
    disconnect(annTableObj, SIGNAL(si_onAnnotationsAdded(QList<Annotation*>)), this, SLOT(sl_lockedStateChanged()));
    annTableList.removeAll(annTableObj);
}

U2SequenceObject* SequenceUndoRedoFramework::getSequenceObject() {
    return dynamic_cast<U2SequenceObject*>(obj);
}

void SequenceUndoRedoFramework::updateObject(QHash<QString, QString> & metaInfo, MaModificationType /*type*/) {// TODO_SVEDIT: why type is not in use?
    U2SequenceObject* seq = getSequenceObject();
    seq->sl_resetDataCaches();
    foreach (AnnotationTableObject* table, annTableList) {
        table->emit_update();
    }
    emit si_updateRequired();

    if (metaInfo.contains("sequence_region")) {
        U2Region region = U2Region::fromString(metaInfo.value("sequence_region"));
        emit si_affectedRegion(region);
    }
}

AnnotationUndoRedoFramework::AnnotationUndoRedoFramework(QObject *p, AnnotationTableObject *annTableObj)
    : UndoRedoFramework(p, annTableObj) {
    connect(annTableObj, SIGNAL(si_onAnnotationModified(AnnotationModification)), SLOT(sl_lockedStateChanged()));
    connect(annTableObj, SIGNAL(si_modifiedStateChanged()), SLOT(sl_lockedStateChanged()));
    connect(annTableObj, SIGNAL(si_nameChanged(QString)), SLOT(sl_lockedStateChanged()));
    connect(annTableObj, SIGNAL(si_onAnnotationsRemoved(QList<Annotation*>)), SLOT(sl_lockedStateChanged()));
    // TODO_SVEDIT: more connection should be added for any kind of editing of annotation

    // TODO_SVEDIT: tmps
    connect(undoAction, SIGNAL(triggered(bool)), SLOT(sl_lockedStateChanged()));
    connect(redoAction, SIGNAL(triggered(bool)), SLOT(sl_lockedStateChanged()));
}

AnnotationTableObject* AnnotationUndoRedoFramework::getAnnotationTableObject() {
    return dynamic_cast<AnnotationTableObject*>(obj);
}

void AnnotationUndoRedoFramework::updateObject(QHash<QString, QString> & , MaModificationType /*type*/) {
    // TODO_SVEDIT: use type??
    AnnotationTableObject* tableObj = getAnnotationTableObject();
    tableObj->emit_update();
}

MsaUndoRedoFramework::MsaUndoRedoFramework(QObject *p, MultipleAlignmentObject *maObj)
    : UndoRedoFramework(p, maObj) {
    connect(maObj, SIGNAL(si_alignmentChanged(const MultipleAlignment&, const MaModificationInfo&)),
                   SLOT(sl_alignmentChanged()));
    connect(maObj, SIGNAL(si_completeStateChanged(bool)), SLOT(sl_completeStateChanged(bool)));
}

void MsaUndoRedoFramework::sl_alignmentChanged() {
    checkUndoRedoEnabled();
}

void MsaUndoRedoFramework::updateObject(QHash<QString, QString> & , MaModificationType type) {
    MaModificationInfo modInfo;
    modInfo.type = type;
    getMaObject()->updateCachedMultipleAlignment(modInfo);
}

MultipleAlignmentObject* MsaUndoRedoFramework::getMaObject() const {
    // TODO_SVEDIT: add safe_point for casting
    return dynamic_cast<MultipleAlignmentObject*>(obj);
}


} // namespace
