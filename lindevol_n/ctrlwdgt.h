#ifndef CTRLWDGT_H
#define CTRLWDGT_H

#include <qwidget.h>
#include <qpushbt.h>

#include "lndglobl.h"


class ControlWidget : public QWidget
{
  Q_OBJECT

public:
  ControlWidget(QWidget *parent = 0, const char *name = 0);
  ~ControlWidget();

public slots:
  void set_finish_flag(int fflag = 1);
  void set_finish_flag();

private:
  QPushButton *quit;
};

#endif

