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

#ifndef _U2_HMMER_BUILD_TASK_H_
#define _U2_HMMER_BUILD_TASK_H_

#include <U2Core/ExternalToolRunTask.h>

#include "HmmerBuildSettings.h"

namespace U2 {

class SaveAlignmentTask;

class HmmerBuildTask : public ExternalToolRunTask {
    Q_OBJECT
public:
    HmmerBuildTask(const HmmerBuildSettings &settings, const QString &stockholmMsaUrl);

    const QString & getHmmProfileUrl() const;
    static QString getReport(const Task *task, const HmmerBuildSettings &settings, const QString &msaUrl);

private:
    void prepare();
    QString generateReport() const;

    static QStringList getArguments(const HmmerBuildSettings &settings, const QString &stockholmMsaUrl);

    HmmerBuildSettings settings;
    const QString stockholmMsaUrl;
};

}   // namespace U2

#endif // _U2_HMMER_BUILD_TASK_H_
