/****************************************************************************
** ControlWidget meta object code from reading C++ file 'ctrlwdgt.h'
**
** Created: Tue Dec 14 18:52:00 1999
**      by: The Qt Meta Object Compiler ($Revision: 1.1.1.1 $)
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

void ControlWidget::initMetaObject()
{
    if ( metaObj )
	return;
    if ( strcmp(QWidget::className(), "QWidget") != 0 )
	badSuperclassWarning("ControlWidget","QWidget");
    if ( !QWidget::metaObject() )
	QWidget::initMetaObject();
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
