/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2019 UniPro <ugene@unipro.ru>
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

#ifndef _U2_GT_UTILS_DASHBOARD_H_
#define _U2_GT_UTILS_DASHBOARD_H_

#include "GTGlobals.h"
#include <primitives/GTWebView.h>

class QWebView;
class QWebElement;
class QTabWidget;

namespace U2 {

class GTUtilsDashboard {
public:
    enum Tabs{Overview, Input, ExternalTools};
    static QWebView* getDashboard(HI::GUITestOpStatus &os);
    static QTabWidget* getTabWidget(HI::GUITestOpStatus &os);

    static QStringList getOutputFiles(HI::GUITestOpStatus &os, const QString &producerName = "");
    static void clickOutputFile(HI::GUITestOpStatus &os, const QString &outputFileName, const QString &producerName = "");

    static HI::HIWebElement findElement(HI::GUITestOpStatus &os, QString text, QString tag = "*", bool exactMatch = false);
    static HI::HIWebElement findTreeElement(HI::GUITestOpStatus &os, QString text);
    static HI::HIWebElement findContextMenuElement(HI::GUITestOpStatus &os, QString text);
    static void click(HI::GUITestOpStatus &os, HI::HIWebElement el, Qt::MouseButton button = Qt::LeftButton);
    static bool areThereNotifications(HI::GUITestOpStatus &os);
    static void openTab(HI::GUITestOpStatus &os, Tabs tab);

private:
    static QMap<QString, Tabs> initTabMap();
    static const QMap<QString, Tabs> tabMap;
};

}   // namespace U2

#endif // _U2_GT_UTILS_DASHBOARD_H_
