#include <QKeyEvent>
#include <QPainter>
#include <QPen>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "hermite.h"
#include "spline.h"
#include "window.h"

#define DEFAULT_A -1.0
#define DEFAULT_B 1.0
#define DEFAULT_N 10
#define DEFAULT_K 0

#define MIN_N 4
#define MAX_K 6
#define SAMPLES 1200

Window::Window(QWidget *parent) : QWidget(parent)
{
    a = DEFAULT_A;
    b = DEFAULT_B;
    n = DEFAULT_N;
    k = DEFAULT_K;

    mode = 0;
    scale_power = 0;
    perturbation = 0;

    x_nodes = NULL;
    f_values = NULL;
    df_values = NULL;
    hermite_coef = NULL;
    spline_coef = NULL;
}

Window::~Window()
{
    clear_data();
}

void Window::clear_data()
{
    delete[] x_nodes;
    delete[] f_values;
    delete[] df_values;
    delete[] hermite_coef;
    delete[] spline_coef;

    x_nodes = NULL;
    f_values = NULL;
    df_values = NULL;
    hermite_coef = NULL;
    spline_coef = NULL;
}

QSize Window::minimumSizeHint() const
{
    return QSize(400, 300);
}

QSize Window::sizeHint() const
{
    return QSize(1000, 700);
}

int Window::parse_command_line(int argc, char *argv[])
{
    if (argc != 5)
        return -1;

    if (sscanf(argv[1], "%lf", &a) != 1)
        return -2;

    if (sscanf(argv[2], "%lf", &b) != 1)
        return -3;

    if (sscanf(argv[3], "%d", &n) != 1)
        return -4;

    if (sscanf(argv[4], "%d", &k) != 1)
        return -5;

    if (b <= a || n < MIN_N || k < 0 || k > MAX_K)
        return -6;

    return rebuild();
}

double Window::function_value(double x) const
{
    switch (k) {
    case 0:
        return 1.0;
    case 1:
        return x;
    case 2:
        return x * x;
    case 3:
        return x * x * x;
    case 4:
        return x * x * x * x;
    case 5:
        return exp(x);
    case 6:
        return 1.0 / (25.0 * x * x + 1.0);
    }

    return 0.0;
}

double Window::derivative_value(double x) const
{
    switch (k) {
    case 0:
        return 0.0;
    case 1:
        return 1.0;
    case 2:
        return 2.0 * x;
    case 3:
        return 3.0 * x * x;
    case 4:
        return 4.0 * x * x * x;
    case 5:
        return exp(x);
    case 6:
        return -50.0 * x / pow(25.0 * x * x + 1.0, 2.0);
    }

    return 0.0;
}

double Window::max_abs_function() const
{
    double max_abs = 0.0;

    for (int i = 0; i < SAMPLES; i++) {
        double x = a + (b - a) * i / (SAMPLES - 1.0);
        double y = fabs(function_value(x));

        if (y > max_abs)
            max_abs = y;
    }

    return max_abs;
}

double Window::max_hermite_error() const
{
    double max_err = 0.0;

    for (int i = 0; i < SAMPLES; i++) {
        double x = a + (b - a) * i / (SAMPLES - 1.0);
        double err = fabs(hermite_compute(x, a, b, n,
                                          x_nodes, hermite_coef) -
                          function_value(x));

        if (err > max_err)
            max_err = err;
    }

    return max_err;
}

double Window::max_spline_error() const
{
    double max_err = 0.0;

    for (int i = 0; i < SAMPLES; i++) {
        double x = a + (b - a) * i / (SAMPLES - 1.0);
        double err = fabs(spline_compute(x, a, b, n,
                                         x_nodes, spline_coef) -
                          function_value(x));

        if (err > max_err)
            max_err = err;
    }

    return max_err;
}

const char *Window::function_name() const
{
    switch (k) {
    case 0:
        return "f(x) = 1";
    case 1:
        return "f(x) = x";
    case 2:
        return "f(x) = x^2";
    case 3:
        return "f(x) = x^3";
    case 4:
        return "f(x) = x^4";
    case 5:
        return "f(x) = exp(x)";
    case 6:
        return "f(x) = 1 / (25x^2 + 1)";
    }

    return "";
}

const char *Window::mode_name() const
{
    switch (mode) {
    case 0:
        return "f и Эрмит";
    case 1:
        return "f и сплайн";
    case 2:
        return "f, Эрмит и сплайн";
    case 3:
        return "погрешности";
    }

    return "";
}

void Window::build_nodes()
{
    double max_abs = max_abs_function();

    for (int i = 0; i < n; i++) {
        x_nodes[i] = a + (b - a) * i / (n - 1.0);
        f_values[i] = function_value(x_nodes[i]);
        df_values[i] = derivative_value(x_nodes[i]);
    }

    f_values[n / 2] += perturbation * 0.1 * max_abs;
}

int Window::rebuild()
{
    clear_data();

    x_nodes = new double[n];
    f_values = new double[n];
    df_values = new double[n];
    hermite_coef = new double[4 * (n - 1)];
    spline_coef = new double[4 * (n - 1)];

    if (x_nodes == NULL || f_values == NULL || df_values == NULL ||
        hermite_coef == NULL || spline_coef == NULL) {
        clear_data();
        return -1;
    }

    build_nodes();

    if (hermite_init(n, x_nodes, f_values, df_values, hermite_coef) != 0)
        return -2;

    if (spline_init(n, x_nodes, f_values, spline_coef) != 0)
        return -3;

    update();

    return 0;
}

double Window::current_left() const
{
    double center = 0.5 * (a + b);
    double half = 0.5 * (b - a) / pow(2.0, scale_power);

    return center - half;
}

double Window::current_right() const
{
    double center = 0.5 * (a + b);
    double half = 0.5 * (b - a) / pow(2.0, scale_power);

    return center + half;
}

double Window::graph_value(double x, int graph_id) const
{
    double f = function_value(x);
    double h = hermite_compute(x, a, b, n, x_nodes, hermite_coef);
    double s = spline_compute(x, a, b, n, x_nodes, spline_coef);

    switch (graph_id) {
    case 0:
        return f;
    case 1:
        return h;
    case 2:
        return s;
    case 3:
        return h - f;
    case 4:
        return s - f;
    }

    return 0.0;
}

void Window::draw_graph(QPainter &painter, double left, double right,
                        int graph_id) const
{
    double x1 = left;
    double y1 = graph_value(x1, graph_id);

    for (int i = 1; i < SAMPLES; i++) {
        double x2 = left + (right - left) * i / (SAMPLES - 1.0);
        double y2 = graph_value(x2, graph_id);

        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));

        x1 = x2;
        y1 = y2;
    }
}

void Window::draw_info(QPainter &painter, double max_abs) const
{
    painter.setPen(Qt::black);

    painter.drawText(10, 20,
                     QString("Функция k = %1, %2")
                         .arg(k)
                         .arg(function_name()));

    painter.drawText(10, 40,
                     QString("Число узлов n = %1, режим: %2")
                         .arg(n)
                         .arg(mode_name()));

    painter.drawText(10, 60,
                     QString("Масштаб s = %1, возмущение p = %2")
                         .arg(scale_power)
                         .arg(perturbation));

    painter.drawText(10, 80, QString("max |F| = %1").arg(max_abs));

    if (mode == 3) {
        painter.drawText(10, 100,
                         QString("max |H-f| = %1")
                             .arg(max_hermite_error()));

        painter.drawText(10, 120,
                         QString("max |S-f| = %1")
                             .arg(max_spline_error()));
    }
}

void Window::paintEvent(QPaintEvent * /* event */)
{
    QPainter painter(this);
    QPen pen;

    double left = current_left();
    double right = current_right();

    double ymin = 0.0;
    double ymax = 0.0;
    int first = 1;

    for (int i = 0; i < SAMPLES; i++) {
        double x = left + (right - left) * i / (SAMPLES - 1.0);

        int graph_ids[3];
        int count = 0;

        if (mode == 0) {
            graph_ids[0] = 0;
            graph_ids[1] = 1;
            count = 2;
        } else if (mode == 1) {
            graph_ids[0] = 0;
            graph_ids[1] = 2;
            count = 2;
        } else if (mode == 2) {
            graph_ids[0] = 0;
            graph_ids[1] = 1;
            graph_ids[2] = 2;
            count = 3;
        } else {
            graph_ids[0] = 3;
            graph_ids[1] = 4;
            count = 2;
        }

        for (int j = 0; j < count; j++) {
            double y = graph_value(x, graph_ids[j]);

            if (first) {
                ymin = y;
                ymax = y;
                first = 0;
            } else {
                if (y < ymin)
                    ymin = y;

                if (y > ymax)
                    ymax = y;
            }
        }
    }

    if (fabs(ymax - ymin) < 1.0e-14) {
        ymin -= 1.0;
        ymax += 1.0;
    }

    double delta = 0.05 * (ymax - ymin);
    ymin -= delta;
    ymax += delta;

    double max_abs = fabs(ymin);

    if (fabs(ymax) > max_abs)
        max_abs = fabs(ymax);

    printf("max |F| = %.16e\n", max_abs);

    painter.save();

    painter.translate(0.5 * width(), 0.5 * height());
    painter.scale(width() / (right - left), -height() / (ymax - ymin));
    painter.translate(-0.5 * (left + right), -0.5 * (ymin + ymax));

    pen.setWidth(0);
    pen.setColor(Qt::black);
    painter.setPen(pen);

    painter.drawLine(QPointF(left, 0.0), QPointF(right, 0.0));
    painter.drawLine(QPointF(0.0, ymin), QPointF(0.0, ymax));

    if (mode == 0 || mode == 2) {
        pen.setColor(Qt::blue);
        painter.setPen(pen);
        draw_graph(painter, left, right, 0);

        pen.setColor(Qt::red);
        painter.setPen(pen);
        draw_graph(painter, left, right, 1);
    }

    if (mode == 1) {
        pen.setColor(Qt::blue);
        painter.setPen(pen);
        draw_graph(painter, left, right, 0);

        pen.setColor(Qt::darkGreen);
        painter.setPen(pen);
        draw_graph(painter, left, right, 2);
    }

    if (mode == 2) {
        pen.setColor(Qt::darkGreen);
        painter.setPen(pen);
        draw_graph(painter, left, right, 2);
    }

    if (mode == 3) {
        pen.setColor(Qt::red);
        painter.setPen(pen);
        draw_graph(painter, left, right, 3);

        pen.setColor(Qt::darkGreen);
        painter.setPen(pen);
        draw_graph(painter, left, right, 4);
    }

    painter.restore();

    draw_info(painter, max_abs);
}

void Window::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_0:
        change_func();
        break;

    case Qt::Key_1:
        change_mode();
        break;

    case Qt::Key_2:
        scale_x_up();
        break;

    case Qt::Key_3:
        scale_x_down();
        break;

    case Qt::Key_4:
        increase_n();
        break;

    case Qt::Key_5:
        decrease_n();
        break;

    case Qt::Key_6:
        increase_perturbation();
        break;

    case Qt::Key_7:
        decrease_perturbation();
        break;

    default:
        QWidget::keyPressEvent(event);
        break;
    }
}

void Window::change_func()
{
    k = (k + 1) % (MAX_K + 1);
    perturbation = 0;
    rebuild();
}

void Window::change_mode()
{
    mode = (mode + 1) % 4;
    update();
}

void Window::scale_x_up()
{
    scale_power++;
    update();
}

void Window::scale_x_down()
{
    scale_power--;
    update();
}

void Window::increase_n()
{
    n *= 2;
    rebuild();
}

void Window::decrease_n()
{
    if (n / 2 >= MIN_N)
        n /= 2;

    rebuild();
}

void Window::increase_perturbation()
{
    perturbation++;
    rebuild();
}

void Window::decrease_perturbation()
{
    perturbation--;
    rebuild();
}
