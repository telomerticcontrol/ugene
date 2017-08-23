# include (U2Gui.pri)

MODULE_ID=U2Gui
include( ../../ugene_lib_common.pri )

UGENE_RELATIVE_DESTDIR = ''

QT += network xml svg sql widgets printsupport
DEFINES+= QT_FATAL_ASSERT BUILDING_U2GUI_DLL
INCLUDEPATH += ../U2Private/src

minQtVersion(5, 4, 0){
    QT += websockets webchannel
}

LIBS += -L../../_release -lU2Core -lU2Formats -lU2Private

!debug_and_release|build_pass {

    CONFIG(debug, debug|release) {
        DESTDIR=../../_debug
        LIBS -= -L../../_release -lU2Core -lU2Formats -lU2Private
        LIBS += -L../../_debug -lU2Cored -lU2Formatsd -lU2Privated
    }

    CONFIG(release, debug|release) {
        DESTDIR=../../_release
    }
}

unix {
    target.path = $$UGENE_INSTALL_DIR/$$UGENE_RELATIVE_DESTDIR
    INSTALLS += target
}

