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

// TODO: rename
#ifndef SANGERFROMACEMAKERTASK_H
#define SANGERFROMACEMAKERTASK_H

#include <U2Core/DNAChromatogram.h>
#include <U2Core/GUrl.h>
#include <U2Core/Task.h>


namespace U2 {

class U2FORMATS_EXPORT SangerFromAceMakerTask : public Task {
    Q_OBJECT
public:
    SangerFromAceMakerTask(const GUrl& aceFile, const QMap<QString, DNAChromatogram>& chromMap,
                           const GUrl& outputFile);
    void run();
private:
    QMap<QString, DNAChromatogram> chromMap;
    GUrl aceFile;
    GUrl outputFile;
};

} // namespace

#endif // SANGERFROMACEMAKERTASK_H
