#include "ctrlwdgt.h"


ControlWidget::ControlWidget(QWidget *parent, const char *name)
{
  quit = new QPushButton("Quit", this, "quit");
  quit->resize(100, 20);
  connect(quit, SIGNAL(clicked()), this, SLOT(set_finish_flag()));
}


ControlWidget::~ControlWidget()
{
}


void ControlWidget::set_finish_flag(int fflag)
{
  if (fflag)
    finish_flag = 1;
}


void ControlWidget::set_finish_flag()
{
  set_finish_flag(1);
}

