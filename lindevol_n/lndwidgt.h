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

