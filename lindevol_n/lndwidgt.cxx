#include <ctype.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <qapp.h>
#include <qobject.h>
#include <qscrbar.h>
#include <qlayout.h>
#include <qpainter.h>

#include "lndglobl.h"
#include "lndwidgt.h"


LndWorldWidget::LndWorldWidget(QWidget *parent, const char *name)
    : QWidget(parent, name)
{
  x_offset = 0;
  y_offset = 0;
  cell_size = 4;
  background_color = new QColor(0, 0, 0);
  cell_color[0] = new QColor(255, 0, 0);
  cell_color[1] = new QColor(0, 255, 0);
  cell_color[2] = new QColor(0, 0, 255);
  cell_color[3] = new QColor(255, 255, 0);
  cell_color[4] = new QColor(255, 0, 255);
  cell_color[5] = new QColor(0, 255, 255);
  setMaximumSize(world_width * cell_size, world_height * cell_size);
  resize(world_width * cell_size, world_height * cell_size);
  setBackgroundColor(*background_color);
}


LndWorldWidget::~LndWorldWidget()
{
  int i;

  delete background_color;
  for (i = 0; i < 6; i++)
    delete cell_color[i];
}


int LndWorldWidget::max_width()
{
  return (world_width * cell_size);
}


int LndWorldWidget::max_height()
{
  return (world_height * cell_size);
}


void LndWorldWidget::show_cell(int xpos, int ypos)
{
  int x, y;
  QPainter p;
  QBrush brush;
  QPen pen;

  x = xpos * cell_size - x_offset;
  y = (world_height - ypos - 1) * cell_size - y_offset;
  p.begin(this);
  if (world[xpos][ypos].plant_no > -1)
  {
    brush.setColor(*cell_color[world[xpos][ypos].plant_no % 6]);
    pen.setColor(*cell_color[world[xpos][ypos].plant_no % 6]);
  }
  else
  {
    brush.setColor(*background_color);
    pen.setColor(*background_color);
  }
  brush.setStyle(SolidPattern);
  p.setPen(pen);
  p.setBrush(brush);
  p.fillRect(x, y, cell_size, cell_size, brush);
  p.end();
}


void LndWorldWidget::paintEvent(QPaintEvent *)
{
  int plant_no, cell_no, x, y;
  QPainter p;
  QBrush brush;
  QPen pen;

  p.begin(this);
  for (plant_no = 0; plant_no < world_width; plant_no++)
  {
    if (plant[plant_no])
    {
      for (cell_no = 0; cell_no < plant[plant_no]->num_cells; cell_no++)
      {
	x = plant[plant_no]->cell[cell_no].x * cell_size - x_offset;
	y = (world_height - plant[plant_no]->cell[cell_no].y - 1) * cell_size - y_offset;
	brush.setColor(*cell_color[plant_no % 6]);
	brush.setStyle(SolidPattern);
	pen.setColor(*cell_color[plant_no % 6]);
	p.setPen(pen);
	p.setBrush(brush);
	p.fillRect(x, y, cell_size, cell_size, brush);
      }
    }
  }
  p.end();
}


void LndWorldWidget::set_xoffset(int x)
{
  x_offset = x;
  update();
}


void LndWorldWidget::set_yoffset(int y)
{
  y_offset = y;
  update();
}


LndWorldScrollWidget::LndWorldScrollWidget(QWidget *parent, const char *name)
    : QWidget(parent, name)
{
  world_widget = new LndWorldWidget(this);
  h_scroll = new QScrollBar(0, world_width - world_widget->width(), 1, world_widget->width(), 0, QScrollBar::Horizontal, this, "h_scroll");
  h_scroll->setFixedHeight(16);
  connect(h_scroll, SIGNAL(valueChanged(int)), world_widget, SLOT(set_xoffset(int)));
  v_scroll = new QScrollBar(0, world_height - world_widget->height(), 1, world_widget->height(), 0, QScrollBar::Vertical, this, "v_scroll");
  v_scroll->setFixedWidth(16);
  connect(v_scroll, SIGNAL(valueChanged(int)), world_widget, SLOT(set_yoffset(int)));
  layout = new QGridLayout(this, 2, 2);
  layout->addWidget(world_widget, 0, 0);
  layout->addWidget(h_scroll, 1, 0);
  layout->addWidget(v_scroll, 0, 1);
  layout->setColStretch(0, 1);
  layout->setRowStretch(0, 1);
  layout->addColSpacing(1, v_scroll->width());
  layout->addRowSpacing(1, h_scroll->height());
  setMaximumSize(world_widget->max_width() + v_scroll->width(), world_widget->max_height() + h_scroll->height());
  resize(world_widget->max_width() + v_scroll->width(), world_widget->max_height() + h_scroll->height());
  /* setGeometry(0, 0, world_widget->width() + v_scroll->width(), world_widget->height() + h_scroll->height()); */
}


LndWorldScrollWidget::~LndWorldScrollWidget()
{
}


void LndWorldScrollWidget::resizeEvent(QResizeEvent *)
{
  h_scroll->setRange(0, world_widget->max_width() - world_widget->width());
  v_scroll->setRange(0, world_widget->max_height() - world_widget->height());
}

