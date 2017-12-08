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

#ifndef _U2_UNDO_REDO_FRAMEWORK_H_
#define _U2_UNDO_REDO_FRAMEWORK_H_

#include <U2Core/MultipleSequenceAlignmentObject.h>

#include <QObject>
#include <QAction>

namespace U2 {

class MultipleAlignmentObject;
class U2SequenceObject;

// TODO_SVEDIT: split headers
class UndoRedoFramework : public QObject {
    Q_OBJECT
public:
    UndoRedoFramework(QObject *p, GObject* obj);

    QAction* getUndoAction() const { return undoAction; }
    QAction* getRedoAction() const { return redoAction; }

private slots:
    void sl_lockedStateChanged();
    void sl_completeStateChanged(bool stateComplete);

    void sl_undo();
    void sl_redo();

protected:
    virtual void updateObject(MaModificationType type) = 0;

protected:
    void checkUndoRedoEnabled();

    GObject*   obj;
    bool       stateComplete;

    QAction*   undoAction;
    QAction*   redoAction;

    qint64     undoStepsAvailable;
    qint64     redoStepsAvailable;
};

class SequenceUndoRedoFramework : public UndoRedoFramework {
    Q_OBJECT
public:
    SequenceUndoRedoFramework(QObject *p, U2SequenceObject* seqObj);

signals:
    void si_updateRequired();

private slots:
    void sl_sequenceChanged();

private:
    U2SequenceObject* getSequenceObject();
    void updateObject(MaModificationType type);
};

class MsaUndoRedoFramework : public UndoRedoFramework {
    Q_OBJECT
public:
    MsaUndoRedoFramework(QObject *p, MultipleAlignmentObject* maObj);

private:
    void updateObject(MaModificationType type);

    MultipleAlignmentObject* getMaObject() const;
private slots:
    void sl_alignmentChanged();
};

} // namespace

#endif
