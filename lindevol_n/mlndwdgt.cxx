/****************************************************************************
** LndWorldWidget meta object code from reading C++ file 'lndwidgt.h'
**
** Created: Tue Dec 14 18:51:57 1999
**      by: The Qt Meta Object Compiler ($Revision: 1.1.1.1 $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#if !defined(Q_MOC_OUTPUT_REVISION)
#define Q_MOC_OUTPUT_REVISION 2
#elif Q_MOC_OUTPUT_REVISION != 2
#error "Moc format conflict - please regenerate all moc files"
#endif

#include "lndwidgt.h"
#include <qmetaobject.h>


const char *LndWorldWidget::className() const
{
    return "LndWorldWidget";
}

QMetaObject *LndWorldWidget::metaObj = 0;

void LndWorldWidget::initMetaObject()
{
    if ( metaObj )
	return;
    if ( strcmp(QWidget::className(), "QWidget") != 0 )
	badSuperclassWarning("LndWorldWidget","QWidget");
    if ( !QWidget::metaObject() )
	QWidget::initMetaObject();
    typedef void(LndWorldWidget::*m1_t0)(int);
    typedef void(LndWorldWidget::*m1_t1)(int);
    m1_t0 v1_0 = &LndWorldWidget::set_xoffset;
    m1_t1 v1_1 = &LndWorldWidget::set_yoffset;
    QMetaData *slot_tbl = new QMetaData[2];
    slot_tbl[0].name = "set_xoffset(int)";
    slot_tbl[1].name = "set_yoffset(int)";
    slot_tbl[0].ptr = *((QMember*)&v1_0);
    slot_tbl[1].ptr = *((QMember*)&v1_1);
    metaObj = new QMetaObject( "LndWorldWidget", "QWidget",
	slot_tbl, 2,
	0, 0 );
}
