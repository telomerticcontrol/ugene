/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2011 UniPro <ugene@unipro.ru>
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

#ifndef CreateExternalProcessDialog_h__
#define CreateExternalProcessDialog_h__

#include "ui/ui_ExternalProcessWorkerDialog.h"
#include <U2Lang/Datatype.h>
#include <U2Lang/Attribute.h>
#include <U2Lang/ConfigurationEditor.h>



namespace U2 {

class ExternalProcessConfig;

class CreateExternalProcessDialog: public QWizard {
    Q_OBJECT
public:
    CreateExternalProcessDialog(QWidget *p = NULL);
    CreateExternalProcessDialog(QWidget *p, ExternalProcessConfig *cfg);
    ExternalProcessConfig* config() const {return cfg;}
    bool validate();

public slots:
    void accept();

private slots:
    void sl_addInput();
    void sl_addOutput();
    void sl_deleteInput();
    void sl_deleteOutput();
    void sl_addAttribute();
    void sl_deleteAttribute();
    void sl_generateTemplateString();
    //void sl_OK();

private:
    Ui::CreateExternalProcessWorkerDialog ui;
    ExternalProcessConfig *cfg;
    bool editing;
    static const int INFO_STRINGS_NUM = 5;
};

}

#endif // CreateExternalProcessDialog_h__
