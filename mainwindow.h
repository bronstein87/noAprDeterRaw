#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>
#include <QFile>
#include <QDir>
#include <QTextStream>
#include <QDebug>
#include <QElapsedTimer>
#include <mathfunc.h>
#include <algorithm>
#include <stdlib.h>
using namespace  BOKZMath;
using namespace std;


struct star
{
    double x;
    double y;
    double bright;
};


struct StarObject
{
    int starNum;
    int objNum;
};

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_fastSectorPushButton_clicked();

    void on_fastThreePushButton_clicked();

private:
    Ui::MainWindow *ui;
    void Pereb();
    void DistCadr();
    void getLoc(const QString& fileName, double focus, quint32 martixSize, double pixelSize);
    void newPereb(StarObject* threeCandidates, int count, QVector <double>& l_st, QVector <double>& m_st, QVector <double>& n_st);

    constexpr const static int  MaxDef = 40;
    void IdentObjects(unsigned short graf[MaxDef][MaxDef], unsigned short Nf, unsigned short* im_save, unsigned short* Nmax);

    constexpr const static int mPM = 48;
    constexpr const static int  NumObj = 10;
    constexpr const static int MinDet = 4;
    constexpr const static double Ro = 206264.806;
    constexpr const static double sec30 = 6875.493164;
    constexpr const static double CosP = 0.98293535;
    constexpr const static int distsCount = 45;
    constexpr const static short groupCount = 70;
    constexpr const static short maxGroupSize = 30;
    int BP[mPM];
    double PM[mPM][mPM];
    double Ncos[NumObj][NumObj];
    double Eps[NumObj][NumObj];
    double EPSILON = 35.0;
    int M_in[NumObj][NumObj];
    double Coord[NumObj][4];
    unsigned short Res[2][MaxDef];
    int M_in_dist[distsCount];
    double Eps_dist[distsCount];
    double NCos_dist[distsCount];
    unsigned short MatCmp[MaxDef][MaxDef];
    float MatDif[MaxDef][MaxDef];
    unsigned short NumIdent[MaxDef];
    unsigned short sum_row[MaxDef], im[MaxDef], im_del[MaxDef];
    //#define sec60 3437.746666

    int kn, ke;    int MaxH = 0;
    int Mnst = 0;
    int ind = 0;
    int Nst = 0;
    unsigned short NumDet = 0;
    int NumDet1 = 0;


};

#endif // MAINWINDOW_H
