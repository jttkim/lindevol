/* $Id: ctrlwdgt.h,v 1.2 2000/01/30 03:10:59 kim Exp $ */
/*
 * $Log: ctrlwdgt.h,v $
 * Revision 1.2  2000/01/30 03:10:59  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

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

