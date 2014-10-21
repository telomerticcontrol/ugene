/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2014 UniPro <ugene@unipro.ru>
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

#include <U2Core/AppContext.h>
#include <U2Core/DNAAlphabet.h>
#include <U2Core/DNATranslation.h>
#include <U2Core/L10n.h>
#include <U2Core/TextUtils.h>
#include <U2Core/U2SafePoints.h>

#include "PrimerStatistics.h"

#include "PrimerGroupBox.h"

namespace U2 {

PrimerGroupBox::PrimerGroupBox(QWidget *parent)
: QWidget(parent)
{
    setupUi(this);
    primerEdit->setValidator(new QRegExpValidator(QRegExp("[acgtACGT]+")));
    reverseComplementButton->setIcon(QIcon(":core/images/todo.png"));
    browseButton->setIcon(QIcon(":core/images/todo.png"));

    connect(primerEdit, SIGNAL(textEdited(const QString &)), SLOT(sl_onPrimerChanged(const QString &)));
    connect(reverseComplementButton, SIGNAL(clicked()), SLOT(sl_translate()));
}

void PrimerGroupBox::setTitle(const QString &title) {
    groupBox->setTitle(title);
}

void PrimerGroupBox::sl_onPrimerChanged(const QString &sequence) {
    QString upSequence = sequence.toUpper();
    primerEdit->setText(upSequence);

    QString characteristics;
    if (!upSequence.isEmpty()) {
        characteristics += getTmString(upSequence) + ", ";
        characteristics += QString::number(upSequence.length()) + tr("-mer");
    }
    characteristicsLabel->setText(characteristics);

    emit si_primerChanged();
}

QString PrimerGroupBox::getSequence() const {
    return primerEdit->text();
}

int PrimerGroupBox::getMismatches() const {
    return mismatchesSpinBox->value();
}

void PrimerGroupBox::sl_translate() {
    const DNAAlphabet *alphabet = AppContext::getDNAAlphabetRegistry()->findById(BaseDNAAlphabetIds::NUCL_DNA_DEFAULT());
    SAFE_POINT(NULL != alphabet, L10N::nullPointerError("DNA Alphabet"), );

    DNATranslation *translator = AppContext::getDNATranslationRegistry()->lookupComplementTranslation(alphabet);
    SAFE_POINT(NULL != translator, L10N::nullPointerError("DNA Translator"), );

    QByteArray sequence = primerEdit->text().toLocal8Bit();
    QByteArray translation(sequence.length(), 0);

    translator->translate(sequence.constData(), sequence.length(), translation.data(), translation.length());
    TextUtils::reverse(translation.data(), translation.length());

    primerEdit->setText(translation);
    sl_onPrimerChanged(translation);
}

QString PrimerGroupBox::getTmString(const QString &sequence) {
    double tm = PrimerStatistics::getMeltingTemperature(sequence.toLocal8Bit());
    QString tmString = QString::number(tm, 'f', 2);
    tmString.remove(QRegExp("\\.?0+$"));
    return tr("Tm = ") + tmString + QString::fromLatin1("\x00B0") + "C";
}

} // U2
