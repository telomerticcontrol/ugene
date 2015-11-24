# include (U2Test.pri)

UGENE_RELATIVE_DESTDIR = ''
MODULE_ID=U2Test
include( ../../ugene_lib_common.pri )

QT += xml gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
DEFINES+= QT_FATAL_ASSERT BUILDING_U2TEST_DLL
LIBS += -L../../_release -lU2Core -lhumimit
INCLUDEPATH += ../../libs_3rdparty/humimit/src

!debug_and_release|build_pass {

    CONFIG(debug, debug|release) {
        DESTDIR=../../_debug
        LIBS -= -L../../_release -lU2Core -lhumimit
        LIBS += -L../../_debug -lU2Cored -lhumimitd
    }

    CONFIG(release, debug|release) {
        DESTDIR=../../_release
    }
}

unix {
    target.path = $$UGENE_INSTALL_DIR/$$UGENE_RELATIVE_DESTDIR
    INSTALLS += target
}

