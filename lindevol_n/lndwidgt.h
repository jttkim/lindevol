/* $Id: lndwidgt.h,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: lndwidgt.h,v $
 * Revision 1.2  2000/01/30 03:11:00  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

#ifndef LNDWIDGT_H
#define LNDWIDGT_H

#include <qwidget.h>
#include <qscrbar.h>
#include <qlayout.h>


class LndWorldWidget : public QWidget
{
  Q_OBJECT
public:
  LndWorldWidget(QWidget *parent = 0, const char *name = 0);
  ~LndWorldWidget();
  void show_cell(int xpos, int ypos);
  int max_width();
  int max_height();

public slots:
  void set_xoffset(int x);
  void set_yoffset(int y);

protected:
  void paintEvent(QPaintEvent *);

private:
  QColor *background_color, *cell_color[6];
  int x_offset, y_offset, cell_size;
};


class LndWorldScrollWidget : public QWidget
{
public:
  LndWorldScrollWidget(QWidget *parent = 0, const char *name = 0);
  ~LndWorldScrollWidget();
  LndWorldWidget *world_widget;

protected:
  void resizeEvent(QResizeEvent *);

private:
  QScrollBar *h_scroll, *v_scroll;
  QGridLayout *layout;
};

#endif

