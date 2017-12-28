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

#include <QFileInfo>

#include <U2Lang/Configuration.h>

#include "DatabaseValidator.h"
#include "KrakenClassifyPrompter.h"
#include "KrakenClassifyWorkerFactory.h"

namespace U2 {
namespace Workflow {

bool DatabaseValidator::validate(const Actor *actor, ProblemList &problemList, const QMap<QString, QString> &) const {
    const QString databaseUrl = actor->getParameter(LocalWorkflow::KrakenClassifyWorkerFactory::DATABASE_ATTR_ID)->getAttributeValueWithoutScript<QString>();
    const bool doesDatabaseDirExist = QFileInfo(databaseUrl).exists();
    CHECK_EXT(doesDatabaseDirExist,
              problemList.append(Problem(LocalWorkflow::KrakenClassifyPrompter::tr("The database folder doesn't exist"), actor->getId())),
              false);

    const QStringList files = QStringList() << "database.kdb"
                                            << "database.idx"
                                            << "taxonomy/nodes.dmp"
                                            << "taxonomy/names.dmp";
    QStringList missedFiles;
    foreach (const QString &file, files) {
        if (!QFileInfo(databaseUrl + "/" + file).exists()) {
            missedFiles << file;
        }
    }

    foreach (const QString &missedFile, missedFiles) {
        problemList.append(Problem(LocalWorkflow::KrakenClassifyPrompter::tr("The mandatory database file doesn't exist: %1").arg(databaseUrl + "/" + missedFile), actor->getId()));
    }
    CHECK(missedFiles.isEmpty(), false);

    return true;
}

}   // namesapce Workflow
}   // namespace U2