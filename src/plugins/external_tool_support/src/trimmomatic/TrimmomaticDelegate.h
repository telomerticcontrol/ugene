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

#ifndef _U2_TRIMMOMATIC_DELEGATE_H_
#define _U2_TRIMMOMATIC_DELEGATE_H_

#include <QStringListModel>

#include <U2Designer/DelegateEditors.h>

#include <U2Gui/SaveDocumentController.h>

#include "TrimmomaticSettingsWidgets.h"
#include "TrimmomaticStepsFactory.h"

#include "ui_TrimmomaticPropertyDialog.h"

namespace U2 {
namespace LocalWorkflow {

class TrimmomaticDelegate : public PropertyDelegate {
    Q_OBJECT
public:
    TrimmomaticDelegate(QObject *parent = 0);

    QVariant getDisplayValue(const QVariant &value) const;
    PropertyDelegate* clone();
    QWidget *createEditor(QWidget *parent, 
                          const QStyleOptionViewItem &option, 
                          const QModelIndex &index) const;
    PropertyWidget *createWizardWidget(U2OpStatus &os, 
                                                  QWidget *parent) const;

    void setEditorData(QWidget *editor, const QModelIndex &index) const;
    void setModelData(QWidget *editor, QAbstractItemModel *model, 
                                       const QModelIndex &index) const;

private slots:
    void sl_commit();
};

class TrimmomaticPropertyWidget : public PropertyWidget {
    Q_OBJECT
public:
    TrimmomaticPropertyWidget(QWidget* parent = NULL, DelegateTags* tags = NULL);

    QVariant value();

public slots:
    void setValue(const QVariant& value);

private slots:
    void sl_showDialog();


private:
    QLineEdit *lineEdit;
    QToolButton *toolButton;
    QString text;
};

class TrimmomaticPropertyDialog : public QDialog, private Ui_TrimmomaticPropertyDialog {
    Q_OBJECT
public:
    TrimmomaticPropertyDialog(const QString &value, QWidget *parent);

    QString getValue() const;

public slots:
    void sl_checkOkEnabled();

private slots:
    void sl_selectionChanged();
    void sl_addStep(QAction* a);
    void sl_moveStepUp();
    void sl_moveStepDown();
    void sl_removeStep();

private:
    void emptySelection();
    void disconnectSelectionChanged();
    void connectSelectionChanged();
    void enableButtons(bool setEnabled);
    static QString defaultDir();

    QList<TrimmomaticBaseController*> steps;
    QWidget* currentWidget;
    TrimmomaticDefaultSettingsWidget* defaultWidget;
    QMenu *menu;
};

}// namespace LocalWorkflow
}// namespace U2

#endif // _U2_TRIMMOMATIC_DELEGATE_H_
