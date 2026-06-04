#include <QAction>
#include <QApplication>
#include <QDebug>
#include <QMainWindow>
#include <QMenuBar>

#include "window.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QMainWindow *window = new QMainWindow;
    QMenuBar *tool_bar = new QMenuBar(window);
    Window *graph_area = new Window(window);
    QAction *action;

    if (graph_area->parse_command_line(argc, argv)) {
        qWarning("Wrong input arguments!");
        return -1;
    }

    action = tool_bar->addAction("Функция: 0", graph_area,
                                 SLOT(change_func()));
    action->setShortcut(QString("0"));

    action = tool_bar->addAction("Режим: 1", graph_area,
                                 SLOT(change_mode()));
    action->setShortcut(QString("1"));

    action = tool_bar->addAction("Увеличить масштаб: 2", graph_area,
                                 SLOT(scale_x_up()));
    action->setShortcut(QString("2"));

    action = tool_bar->addAction("Уменьшить масштаб: 3", graph_area,
                                 SLOT(scale_x_down()));
    action->setShortcut(QString("3"));

    action = tool_bar->addAction("Увеличить n: 4", graph_area,
                                 SLOT(increase_n()));
    action->setShortcut(QString("4"));

    action = tool_bar->addAction("Уменьшить n: 5", graph_area,
                                 SLOT(decrease_n()));
    action->setShortcut(QString("5"));

    action = tool_bar->addAction("Увеличить возмущение: 6", graph_area,
                                 SLOT(increase_perturbation()));
    action->setShortcut(QString("6"));

    action = tool_bar->addAction("Уменьшить возмущение: 7", graph_area,
                                 SLOT(decrease_perturbation()));
    action->setShortcut(QString("7"));

    action = tool_bar->addAction("Выход", window, SLOT(close()));
    action->setShortcut(QString("Ctrl+X"));

    tool_bar->setMaximumHeight(30);

    window->setMenuBar(tool_bar);
    window->setCentralWidget(graph_area);
    window->setWindowTitle("Приближение функций");

    window->show();

    return app.exec();
}
