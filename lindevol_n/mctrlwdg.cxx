/****************************************************************************
** ControlWidget meta object code from reading C++ file 'ctrlwdgt.h'
**
** Created: Sun Jan 30 03:21:09 2000
**      by: The Qt Meta Object Compiler ($Revision: 1.2 $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#if !defined(Q_MOC_OUTPUT_REVISION)
#define Q_MOC_OUTPUT_REVISION 2
#elif Q_MOC_OUTPUT_REVISION != 2
#error "Moc format conflict - please regenerate all moc files"
#endif

#include "ctrlwdgt.h"
#include <qmetaobject.h>


const char *ControlWidget::className() const
{
    return "ControlWidget";
}

QMetaObject *ControlWidget::metaObj = 0;


#if QT_VERSION >= 200
static QMetaObjectInit init_ControlWidget(&ControlWidget::staticMetaObject);

#endif

void ControlWidget::initMetaObject()
{
    if ( metaObj )
	return;
    if ( strcmp(QWidget::className(), "QWidget") != 0 )
	badSuperclassWarning("ControlWidget","QWidget");

#if QT_VERSION >= 200
    staticMetaObject();
}

void ControlWidget::staticMetaObject()
{
    if ( metaObj )
	return;
    QWidget::staticMetaObject();
#else

    QWidget::initMetaObject();
#endif

    typedef void(ControlWidget::*m1_t0)(int);
    typedef void(ControlWidget::*m1_t1)();
    m1_t0 v1_0 = &ControlWidget::set_finish_flag;
    m1_t1 v1_1 = &ControlWidget::set_finish_flag;
    QMetaData *slot_tbl = new QMetaData[2];
    slot_tbl[0].name = "set_finish_flag(int)";
    slot_tbl[1].name = "set_finish_flag()";
    slot_tbl[0].ptr = *((QMember*)&v1_0);
    slot_tbl[1].ptr = *((QMember*)&v1_1);
    metaObj = new QMetaObject( "ControlWidget", "QWidget",
	slot_tbl, 2,
	0, 0 );
}
