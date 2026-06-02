#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>

class Window : public QWidget
{
    Q_OBJECT

  private:
    double a;
    double b;
    int n;
    int k;
    int mode;
    int scale_power;
    int perturbation;

    double *x_nodes;
    double *f_values;
    double *df_values;
    double *hermite_coef;
    double *spline_coef;

    int rebuild();
    void clear_data();
    void build_nodes();

    double function_value(double x) const;
    double derivative_value(double x) const;
    double max_abs_function() const;

    double current_left() const;
    double current_right() const;

    double graph_value(double x, int graph_id) const;
    void draw_graph(QPainter &painter, double left, double right,
                    int graph_id) const;
    void draw_info(QPainter &painter, double max_abs) const;

    const char *function_name() const;
    const char *mode_name() const;

  public:
    Window(QWidget *parent);
    ~Window();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;

    int parse_command_line(int argc, char *argv[]);

  public slots:
    void change_func();
    void change_mode();
    void scale_x_up();
    void scale_x_down();
    void increase_n();
    void decrease_n();
    void increase_perturbation();
    void decrease_perturbation();

  protected:
    void paintEvent(QPaintEvent *event);
    void keyPressEvent(QKeyEvent *event);
};

#endif
