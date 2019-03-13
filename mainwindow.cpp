#include "mainwindow.h"
#include "ui_mainwindow.h"



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}




void setBit( char A[],  int k )
{
    A[k/8] |= 1 << (k%8);  // Set the bit at the k-th position in A[i]
}

void clearBit( char A[],  int k )
{
    A[k/8] &= ~(1 << (k%8));
}

bool testBit( char A[],  int k )
{
    return ( (A[k/8] & (1 << (k%8) )) != 0 ) ;
}


void MainWindow::getLoc(const QString& fileName, double focus, quint32 martixSize, double pixelSize)
{
    // focus = focus / pixelSize;
    QFile file (fileName);
    if (file.open(QFile::ReadOnly)) {
        QTextStream in(&file);
        QString currentLine;
        QVector <QPointF> starsCoords;
        QVector <star> stars;
        bool coordsReaded = false;
        while  (in.readLineInto(&currentLine)) {
            if (currentLine.contains("Max")) {
                for (quint32 i = 0; i < NumObj; i++) {
                    in.readLineInto(&currentLine);
                    QStringList list = currentLine.split("  ", QString::SkipEmptyParts);
                    star s;
                    //                    if (i < 6)
                    //                    {
                    //                        s.x = ((static_cast <float> (rand()) / static_cast <float> (RAND_MAX / 2048)) - martixSize / 2.) * pixelSize;
                    //                        s.y = ((static_cast <float> (rand()) / static_cast <float> (RAND_MAX / 2048)) - martixSize / 2.) * pixelSize;
                    //                    }
                    //                    else
                    //                    {
                    //qDebug() << ((static_cast <float> (rand()) / static_cast <float> (RAND_MAX / 2048)) - martixSize / 2.);
                    s.x = (list[2].toDouble() - martixSize / 2.) * pixelSize;
                    s.y = (list[3].toDouble() - martixSize / 2.) * pixelSize;

                    //                    }
                    stars.append(s);
                }
                coordsReaded = true;
            }
            if (currentLine.contains("Матрица ориентации:") && coordsReaded) {
                double orientMatrix[3][3];
                for (int i = 0; i < 3; i ++) {
                    for (int j = 0; j < 3; j ++) {
                        in >> orientMatrix[i][j];
                    }
                }
                break;
            }
        }

        if  (!coordsReaded) {
            qDebug() << starsCoords << " - неполный протокол";
        }

        std::sort(stars.begin(), stars.end(), [](auto& a, auto& b){return a.bright > b.bright;});
        for (const auto& i : stars)
        {
            starsCoords.append(QPointF(i.x, i.y));
        }
        // starsCoords.append(QPointF((list[2].toDouble() - martixSize / 2.) * pixelSize,
        //                   (list[3].toDouble() - martixSize / 2.) * pixelSize));

        for (int i = 0; i < NumObj; i++)
        {
            double length = sqrt(starsCoords[i].x() * starsCoords[i].x() + starsCoords[i].y() * starsCoords[i].y() + focus * focus);
            Coord[i][0] = - starsCoords[i].x() / length;
            Coord[i][1] = - starsCoords[i].y() / length;
            Coord[i][2] = focus / length;
            Coord[i][3] = 0;
        }
    }
}


#define Dist(li,mi,ni,lj,mj,nj) li*lj+mi*mj+ni*nj;
void MainWindow::DistCadr()
{
    short i, j, k;
    double an, res;
    short count = 0;
    for (i = 0; i < NumObj - 1; i++)
    {
        for (j = (i+1); j < NumObj; j++)
        {
            Ncos[i][j] = Dist(Coord[i][0], Coord[i][1], Coord[i][2],
                    Coord[j][0], Coord[j][1], Coord[j][2]);
            float res = 1.0 - Ncos[i][j];
            if (res>0)
                res = sqrt(res * 2.0) * EPSILON / Ro;
            else res = 0.0;
            Eps[i][j] = res;

        }
    }

    for (i = 1; i < NumObj; i++)
    {
        for (j = 0; j < i; j++)
        {
            Ncos[i][j] = Ncos[j][i];
            Eps[i][j] = Eps[j][i];
        }
    }

    for (i = 0; i < NumObj; i++)
    {
        for (j = 0; j < NumObj; j++)
        {
            M_in[i][j] = -10;
        }
    }

    for (i = 0; i < distsCount; i++)
    {
        M_in_dist[i] = -10;
    }
    for (i = 0; i<NumObj - 1; i++)
    {
        for (j = i+1; j < NumObj; j++)
        {
            if ((fabs(Ncos[i][j])<1.0) && (Ncos[i][j] > CosP))
            {
                an = acos((double)Ncos[i][j]) * sec30;
                k = (short)(an);
                M_in[i][j] = k;
            }
            //else M_in[i][j]=-10;
            Eps_dist[count] = Eps[i][j];
            NCos_dist[count] = Ncos[i][j];
            M_in_dist[count++] = M_in[i][j];
        } // end for j
    }
}


void MainWindow::Pereb(void)
{
    short i, j, k, Im, Jm, Km, min1, min2, min3;
    short ii, jj, i1, j1;
    float Pq;
    unsigned short OK, OK1, Bu;

    OK=0;
    i = 0; min3 = MinDet - 1; min2 = MinDet - 2; min1 = MinDet - 3;
    while ((i < (NumObj - min3)) && (OK == 0))        // 1 index by objects
    {
        j=i+1;
        while ((j<(NumObj-min2)) && (OK==0))       // 2 index by objects
        {
            Im=0;
            while ((Im<Mnst) && (OK==0))           // 1 index by objects
            {
                Jm=0;
                while ((Jm<Mnst) && (OK==0))         // 2 index by objects
                {
                    Pq=fabs(PM[Im][Jm]-Ncos[i][j]);
                    if ((Pq<Eps[i][j]) && (Im!=Jm))
                    {
                        k=j+1;
                        while ((k<(NumObj-min1)) && (OK==0))   // 3 index by objects
                        {
                            Km=0;
                            while ((Km<Mnst) && (OK==0))   // 3 index by objects
                            {
                                Pq=fabs(PM[Jm][Km]-Ncos[j][k]);
                                if ((Pq<Eps[j][k]) && (Jm!=Km))
                                {
                                    Pq = fabs(PM[Im][Km] - Ncos[i][k]);
                                    if ((Pq < Eps[i][k]) && (Im != Km))
                                    {
                                        Res[0][0] = i; Res[1][0] = Im;
                                        Res[0][1] = j; Res[1][1] = Jm;
                                        Res[0][2] = k; Res[1][2] = Km;
                                        NumDet = 3; jj = k + 1; ii = 3;
                                        while (jj < NumObj)  //
                                        {
                                            j1 = 0; OK1 = 0;
                                            while ((j1 < Mnst) && (OK1 == 0))
                                            {
                                                Bu = 0;
                                                for (i1 = 0; i1 < ii; i1++)
                                                {
                                                    if (Res[1][i1] == j1) Bu = 1;
                                                }
                                                if (Bu==0)
                                                {
                                                    for (i1 = 0; i1 < ii; i1++)
                                                    {
                                                        Pq = fabs(PM[Res[1][i1]][j1] - Ncos[Res[0][i1]][jj]);
                                                        if (Pq >= Eps[Res[0][i1]][jj])
                                                        {
                                                            Bu = 1;
                                                        }
                                                    }
                                                    if (Bu == 0)
                                                    {
                                                        Res[0][ii] = jj; Res[1][ii] = j1;
                                                        ii++; NumDet++; OK1 = 1;
                                                    }
                                                }
                                                j1++;
                                            }
                                            jj++;
                                        }
                                        if (NumDet >= MinDet)
                                        {
                                            OK=1;
                                        }
                                    }
                                }
                                Km++;
                            }
                            k++;
                        }
                    }
                    Jm++;
                }
                Im++;
            }
            j++;
        }
        i++;
    }
}


void MainWindow::on_fastSectorPushButton_clicked()
{
    struct sectors // каталог секторов
    {
        float l;
        float m;
        float n;
        quint16 count_in_sector; // ВНИМАТЕЛЬНО, ПОД ГАЙА 32
        qint16 shift;
    };

    QVector <sectors> sec;
    QVector <qint16> num;

    QVector<double> l_st;
    QVector<double> m_st;
    QVector<double> n_st;
    QVector <double> angle_arccos;

    QFile lFile("C:/Users/Yumatov/Documents/catalog/Выборочные каталоги/МР/5153_L.CAT");
    if (lFile.open(QIODevice::ReadOnly))
    {
        QDataStream in(&lFile);
        float l;
        while (in.readRawData((char*)(&l), sizeof(float)))
        {
            l_st.append(l);
        }
    }
    QFile mFile("C:/Users/Yumatov/Documents/catalog/Выборочные каталоги/МР/5153_M.CAT");
    if (mFile.open(QIODevice::ReadOnly))
    {
        QDataStream in(&mFile);
        float m;
        while (in.readRawData((char*)(&m), sizeof(float)))
        {
            m_st.append(m);
        }
    }
    QFile nFile("C:/Users/Yumatov/Documents/catalog/Выборочные каталоги/МР/5153_N.CAT");
    if (nFile.open(QIODevice::ReadOnly))
    {
        QDataStream in(&nFile);
        float n;
        while (in.readRawData((char*)(&n), sizeof(float)))
        {
            n_st.append(n);
        }
    }

    QFile sFile("C:/Users/Yumatov/Documents/catalog/Выборочные каталоги/МР/5153REAL_SEC.CAT");
    if (sFile.open(QIODevice::ReadOnly))
    {
        QDataStream in(&sFile);
        sectors s;
        while (in.readRawData((char*)(&s), sizeof(sectors)))
        {
            sec.append(s);
        }
    }

    QFile nuFile("C:/Users/Yumatov/Documents/catalog/Выборочные каталоги/МР/5153REAL_NUN.CAT");
    if (nuFile.open(QIODevice::ReadOnly))
    {
        QDataStream in(&nuFile);
        qint16 n;
        while (in.readRawData((char*)(&n), sizeof(short)))
        {
            num.append(n);
        }
    }

    const int intervalCount = 1270; // 635
    const int bitCount = 160;
    const int sCount = 5153;
    double delta = 30.0;
    char histogramm[sCount][bitCount];
    for (int i = 0; i < sCount; i++)
    {
        for (int j = 0; j < bitCount; j++)
        {
            histogramm[i][j] = 0;
        }
    }
    for (int i = 0; i < sec.size(); i++)
    {
        for (int j = sec[i].shift; j < sec[i].shift + sec[i].count_in_sector - 1; j++)
        {
            for (int k = j + 1; k < sec[i].shift + sec[i].count_in_sector; k++)
            {
                double scalarProduct = l_st[num[j]] * l_st[num[k]] + m_st[num[j]] * m_st[num[k]] + n_st[num[j]] * n_st[num[k]];
                double angle = acosm(scalarProduct) * radToDegrees * 3600;
                int pos = angle / delta;

                if (pos < intervalCount)
                {
                    setBit(histogramm[num[j]], pos);
                    setBit(histogramm[num[k]], pos);
                    if (pos != 0)
                    {
                        setBit(histogramm[num[j]], pos - 1);
                        setBit(histogramm[num[k]], pos - 1);
                        if (pos > 1)
                        {
                            setBit(histogramm[num[j]], pos - 2);
                            setBit(histogramm[num[k]], pos - 2);
                        }
                    }
                    if (pos < intervalCount - 1)
                    {
                        setBit(histogramm[num[j]], pos + 1);
                        setBit(histogramm[num[k]], pos + 1);
                        if (pos < intervalCount - 2)
                        {
                            setBit(histogramm[num[j]], pos + 2);
                            setBit(histogramm[num[k]], pos + 2);
                        }
                    }
                }
            }
        }
    }

    QFile histBin("hist.CAT");
    histBin.open(QIODevice::WriteOnly);
    QDataStream out(&histBin);
    for (int i = 0; i < sCount; i++)
    {
        out.writeRawData((char*)histogramm[i], bitCount * sizeof(char));
    }
    QDir testDir("C:/Users/Yumatov/Documents/build-catalog-Desktop_Qt_5_8_0_MinGW_32bit3-Release/test");
    auto files = testDir.entryList(QDir::Files);
    double focus = 59.782900;
    //double focus = 60;

    const int maxSectorSize = 200;
    int hist [NumObj][maxSectorSize];
    for (const auto& f : files)
    {
        NumDet1 = 0;
        getLoc(testDir.path() + "/" + f, focus, 2048, 0.0055);
        QElapsedTimer timer;
        timer.start();
        DistCadr();
        for (int s = 0; s < sec.size(); s++)
        {
            Nst = sec[s].count_in_sector;
            for (int i = 0; i < 10; i++)
            {
                for (int j = 0; j < 200; j++)
                {
                    hist[i][j] = 0;
                }
            }
            for (int i = 0; i<(NumObj - 1); i++)
            {
                for (int j = (i + 1); j < NumObj; j++)
                {
                    for (int si = sec[s].shift; si < sec[s].shift + sec[s].count_in_sector; si++)
                    {
                        if (testBit(histogramm[num[si]], M_in[i][j]))
                        {
                            ++hist[i][si - sec[s].shift];
                            ++hist[j][si - sec[s].shift];
                        }
                    }
                }
            }

            //                        qDebug() << "\nHIST\n";
            //                        for (int i = 0; i < NumObj; i++)
            //                        {
            //                            QString line;
            //                            //sort(hist[i], hist[i] + 200, [](auto& a, auto& b){return a > b;});
            //                            for (int j = 0; j < 200; j++)
            //                            {
            //                                line.append(QString::number(hist[i][j]) + " ");
            //                            }
            //                            line.append("\n");
            //                            qDebug() << line;
            //                        }
            //                        qDebug() << "\nHIST_END";

            Mnst = 0;
            for (int i = 0; i < NumObj; i++)
            {
                MaxH = 0;
                bool exit = false;
                for (int i1 = 0; i1 < Nst; i1++)
                {
                    if (hist[i][i1] > MaxH)
                    {
                        MaxH = hist[i][i1];
                    }
                }
                if (MaxH >= MinDet)
                {
                    int ih = 0;
                    for (int i1 = 0; i1 < Nst; i1++)
                    {
                        if (hist[i][i1] == MaxH)
                        {
                            BP[Mnst + ih] = i1;
                            ih++;
                            if (Mnst + ih >= mPM)
                            {
                                exit = true;
                            }
                        }
                    }
                    if (exit)
                    {
                        break;
                    }
                    Mnst += ih; NumDet++;
                }
            }
            if (NumDet >= MinDet)  //HO/TO MinDet=4
            {
                for (int i = 0; i < Mnst; i++)
                {
                    for (int j1 = 0; j1 < Mnst; j1++)
                    {
                        PM[i][j1] = l_st[num[sec[s].shift + BP[i]]] * l_st[num[sec[s].shift + BP[j1]]]
                                +m_st[num[sec[s].shift + BP[i]]] * m_st[num[sec[s].shift + BP[j1]]]
                                +n_st[num[sec[s].shift + BP[i]]] * n_st[num[sec[s].shift + BP[j1]]];
                    }
                }
                NumDet = 0;
                Pereb();
                if(NumDet >= MinDet)
                {
                    NumDet1 = NumDet;
                }
            }
        }
        qDebug() << timer.nsecsElapsed() << NumDet1;
    }
}

#pragma pack(push,1)
struct RangeAttrs
{
    int start;
    int count;
};
#pragma pack(pop)


struct Candidate
{
    short dist;
    char indexArray;
    char firstIndex;
    char secondIndex;
};


#pragma pack(push,1)
struct StarPair
{
    short fStar;
    short sStar;
};
#pragma pack(pop)



int compareCandidate (const void* p1, const void* p2)
{
    const Candidate* f = (const Candidate*)p1;
    const Candidate* s = (const Candidate*)p2;
    if (f->dist < s->dist) return -1;
    else if (f->dist > s->dist) return 1;
    else return 0;
}

int compareStarObject (const void* p1, const void* p2)
{

    StarObject* tmp1 = (StarObject*)p1;
    StarObject* tmp2 = (StarObject*)p2;

    if ( (tmp1->starNum <  tmp2->starNum)
         || (tmp1->starNum ==  tmp2->starNum
             && (tmp1->objNum <  tmp2->objNum)))
    {
        return -1;
    }
    if ((tmp1->starNum ==  tmp2->starNum
         && (tmp1->objNum ==  tmp2->objNum)))
    {
        return 0;
    }
    if ( (tmp1->starNum >  tmp2->starNum)
         || (tmp1->starNum ==  tmp2->starNum
             && (tmp1->objNum >  tmp2->objNum)))
    {
        return 1;
    }
}


short uniques(StarObject array[], size_t size)
{
    short newSize = 0;
    if (size == 0)
        return 0;

    StarObject* result = array;
    StarObject* first = array;
    StarObject* last = (array + size);
    while (++first != last) {
        if (!(result->starNum == first->starNum)
                || (result->starNum == first->starNum
                    && result->objNum != first->objNum)) {
            *(++result) = *first;
            ++newSize;
        }
    }
    return ++newSize;
}

#if defined(__i386__)

static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long int x;
    __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
    return x;
}

#elif defined(__x86_64__)

static __inline__ unsigned long long rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#endif

void assingObjectIndexToStar(Candidate* mainp, Candidate* addp, char* sObjects)
{
    bool firstFirst = false;
    if ((firstFirst = (mainp->firstIndex == addp->firstIndex))
            || mainp->firstIndex == addp->secondIndex)
    {
        sObjects[0] =  mainp->firstIndex;
        sObjects[1] = mainp->secondIndex;
        sObjects[2] = firstFirst ? addp->secondIndex : addp->firstIndex;
        return;
    }
    sObjects[0] =  mainp->secondIndex;
    sObjects[1] = mainp->firstIndex;
    sObjects[2] = mainp->secondIndex == addp->firstIndex ? addp->secondIndex : addp->firstIndex;
}

unsigned __int64 time1, time2;

void MainWindow::on_fastThreePushButton_clicked()
{
    QVector <StarPair> nCatalog;
    QVector <RangeAttrs> shifts;

    const  short intervalCount = 1270;
    QVector <QVector <StarPair>> rawCatalog(intervalCount);
    QVector<double> l_st;
    QVector<double> m_st;
    QVector<double> n_st;
    QVector <double> angle_arccos;

    QFile lFile("C:/Users/Yumatov/Documents/catalog/Выборочные каталоги/МР/5153_L.CAT");
    if (lFile.open(QIODevice::ReadOnly))
    {
        QDataStream in(&lFile);
        float l;
        while (in.readRawData((char*)(&l), sizeof(float)))
        {
            l_st.append(l);
        }
    }
    QFile mFile("C:/Users/Yumatov/Documents/catalog/Выборочные каталоги/МР/5153_M.CAT");
    if (mFile.open(QIODevice::ReadOnly))
    {
        QDataStream in(&mFile);
        float m;
        while (in.readRawData((char*)(&m), sizeof(float)))
        {
            m_st.append(m);
        }
    }
    QFile nFile("C:/Users/Yumatov/Documents/catalog/Выборочные каталоги/МР/5153_N.CAT");
    if (nFile.open(QIODevice::ReadOnly))
    {
        QDataStream in(&nFile);
        float n;
        while (in.readRawData((char*)(&n), sizeof(float)))
        {
            n_st.append(n);
        }
    }

    double delta = 30.0;

    for (int i = 0; i < l_st.size() - 1; i++)
    {
        for (int j = i + 1; j < l_st.size(); j++)
        {
            double scalarProduct = l_st[i] * l_st[j] + m_st[i] * m_st[j] + n_st[i] * n_st[j];
            double angle = acosm(scalarProduct) * radToDegrees * 3600;
            int pos = angle / delta;
            if (pos < intervalCount)
            {
                StarPair p;
                p.fStar = i;
                p.sStar = j;
                rawCatalog[pos].append(p);
            }
        }
    }
    int shift = 0;
    for (int i = 0; i < intervalCount; i++)
    {
        for (const auto& j : rawCatalog[i])
        {
            nCatalog.append(j);
        }
        RangeAttrs attrs;
        attrs.start = shift;
        attrs.count = rawCatalog[i].size();
        shifts.append(attrs);
        shift += rawCatalog[i].size();
    }

    QFile shiftsFile("SHIFTS.CAT");
    shiftsFile.open(QIODevice::WriteOnly);
    QDataStream out(&shiftsFile);
    for (int i = 0; i < shifts.size(); i++)
    {
        out.writeRawData((char*)&shifts[i], sizeof(RangeAttrs));
    }
    shiftsFile.close();

    QFile histFile("HIST.CAT");
    histFile.open(QIODevice::WriteOnly);
    QDataStream out1(&histFile);
    for (int i = 0; i < nCatalog.size(); i++)
    {
        out1.writeRawData((char*)&nCatalog[i], sizeof(StarPair));
    }
    histFile.close();

    //qDebug() << "qq";
    //QDir testDir("C:/Users/Yumatov/Documents/build-catalog-Desktop_Qt_5_8_0_MinGW_32bit3-Release/test");
    //QDir testDir("C:/Users/Yumatov/Documents/build-catalog-Desktop_Qt_5_8_0_MinGW_32bit3-Release/ttest");
    //QDir testDir("C:/Users/Yumatov/Documents/build-catalog-Desktop_Qt_5_8_0_MinGW_32bit3-Release/tttest");
    QDir testDir("C:/Users/Yumatov/Documents/build-catalog-Desktop_Qt_5_8_0_MinGW_32bit3-Release/ttest");
    //QDir testDir("C:/Users/Yumatov/Documents/build-catalog-Desktop_Qt_5_8_0_MinGW_32bit3-Release/nnew_test");



    auto files = testDir.entryList(QDir::Files);
    double focus = 59.782900;
    //double focus = 60;
    // double focus = 59.31;
    //double MApr [3] = {-0.3842312,   0.6449825 , -0.6605785};

    StarObject tmpStar1;
    StarObject tmpStar2;
    StarObject tmpStar3;
    int cnt = 0;

    char sObjects[3];
    //short* tmpStars[3] {&tmpStar1, &tmpStar2, &tmpStar3};
    QVector <quint64> times;
    int count = 0;
    for (const auto& f : files)
    {
        NumDet = 0;
        qDebug() << f;
        StarObject threeCandidates[groupCount][maxGroupSize];
        short threeCandidatesCount[groupCount];
        char hist [5153][6];
        getLoc(testDir.path() + "/" + f, focus, 2048, 0.0055);
        QElapsedTimer timer;
        //time1 = rdtsc();
        timer.start();
        memset(hist, 0, sizeof(hist[0][0]) * 5153 * 6);
        DistCadr();
        for (int i = 0; i < distsCount; i++)
        {
            if (M_in_dist[i] >= 0)
            {
                for (int j = -2; j < 3; j++)
                {
                    int pos = M_in_dist[i]  - j;
                    if (pos < 0 || pos >= intervalCount)
                    {
                        continue;
                    }
                    for (int k = shifts[pos].start; k < shifts[pos].start + shifts[pos].count; k++)
                    {
                        double cos = l_st[nCatalog[k].fStar] * l_st[nCatalog[k].sStar]
                                + m_st[nCatalog[k].fStar] * m_st[nCatalog[k].sStar]
                                + n_st[nCatalog[k].fStar] * n_st[nCatalog[k].sStar];
                        if (fabs(cos - NCos_dist[i]) <= Eps_dist[i])
                        {
                            setBit(hist[nCatalog[k].fStar], i);
                            setBit(hist[nCatalog[k].sStar], i);
                        }
                    }
                }
            }
        }


        memset(threeCandidatesCount, 0, sizeof(short) * groupCount);
        // time2 = rdtsc();
        // times.append(time2 - time1);
        // time1 = rdtsc();
        //qint32 count = 0;
        qint32 distIndexFirstStar = -1;
        for (int i = 0; i < NumObj - 2; i++)
        {
            short offsetI = 1;
            for (int t = i; t >= 0; t--)
            {
                offsetI += t;
            }
            for (int j = i + 1; j < NumObj - 1; j++)
            {
                short offsetJ = 1;
                for (int t = j; t >= 0; t--)
                {
                    offsetJ += t;
                }
                Candidate candidates[3];
                for (int k = j + 1; k < NumObj; k++)
                {
                    if (M_in[i][j] < 0 || M_in[i][k] < 0 || M_in[j][k] < 0)
                    {
                        continue;
                    }
                    candidates[0].dist = M_in[i][j];
                    candidates[0].indexArray = (i * NumObj) - offsetI + (j - i);
                    candidates[0].firstIndex = i;
                    candidates[0].secondIndex = j;

                    candidates[1].dist = M_in[i][k];
                    candidates[1].indexArray = (i * NumObj) - offsetI + (k - i);
                    candidates[1].firstIndex = i;
                    candidates[1].secondIndex = k;

                    candidates[2].dist = M_in[j][k];
                    candidates[2].indexArray = (j * NumObj) - offsetJ + (k - j);
                    candidates[2].firstIndex = j;
                    candidates[2].secondIndex = k;

                    qsort(candidates, 3, sizeof(Candidate), compareCandidate);
                    short min1 = candidates[0].dist - 2 >= 0 ? candidates[0].dist - 2
                            : candidates[0].dist - 1 >= 0 ? candidates[0].dist - 1 : candidates[0].dist;

                    short max1 = candidates[0].dist + 2 < intervalCount ? candidates[0].dist + 2
                            : candidates[0].dist - 1 < intervalCount ? candidates[0].dist + 1 : candidates[0].dist;

                    for (int d1 = min1; d1 < max1 + 1; d1++)
                    {
                        for (int p = shifts[d1].start; p < shifts[d1].start + shifts[d1].count; p++)
                        {
                            tmpStar1.starNum = -1;
                            tmpStar2.starNum = -1;
                            tmpStar3.starNum = -1;
                            distIndexFirstStar = -1;

                            if (testBit(hist[nCatalog[p].fStar], candidates[1].indexArray)
                                    && testBit(hist[nCatalog[p].sStar], candidates[2].indexArray))
                            {
                                distIndexFirstStar = 1;
                            }
                            else if (testBit(hist[nCatalog[p].fStar], candidates[2].indexArray)
                                     && testBit(hist[nCatalog[p].sStar], candidates[1].indexArray))
                            {
                                distIndexFirstStar = 2;
                            }
                            if (distIndexFirstStar != -1)
                            {
                                tmpStar1.starNum = nCatalog[p].fStar;
                                tmpStar2.starNum = nCatalog[p].sStar;
                            }

                            if (tmpStar1.starNum != -1)
                            {
                                short min2 = candidates[1].dist - 2 >= 0 ? candidates[1].dist - 2
                                        : candidates[1].dist - 1 >= 0 ? candidates[1].dist - 1 : candidates[1].dist;

                                short max2 = candidates[1].dist + 2 < intervalCount ? candidates[1].dist + 2
                                        : candidates[1].dist - 1 < intervalCount ? candidates[1].dist + 1 : candidates[1].dist;

                                bool threeClosed = false;
                                for (int d2 = min2; d2 < max2 + 1; d2++)
                                {
                                    if (threeClosed) break;
                                    for (int pm = shifts[d2].start; pm < shifts[d2].start + shifts[d2].count; pm++)
                                    {
                                        if ((testBit(hist[nCatalog[pm].fStar], candidates[0].indexArray)
                                             && testBit(hist[nCatalog[pm].sStar], candidates[2].indexArray))
                                                || (testBit(hist[nCatalog[pm].fStar], candidates[2].indexArray)
                                                    && testBit(hist[nCatalog[pm].sStar], candidates[0].indexArray)))
                                        {
                                            if ((nCatalog[pm].fStar == tmpStar1.starNum
                                                 && nCatalog[pm].sStar != tmpStar2.starNum)
                                                    || (nCatalog[pm].fStar == tmpStar2.starNum
                                                        && nCatalog[pm].sStar != tmpStar1.starNum))
                                            {
                                                tmpStar3.starNum = nCatalog[pm].sStar;
                                                threeClosed = true;
                                                break;
                                            }
                                            if ((nCatalog[pm].sStar == tmpStar1.starNum
                                                 && nCatalog[pm].fStar != tmpStar2.starNum)
                                                    || (nCatalog[pm].sStar == tmpStar2.starNum
                                                        && nCatalog[pm].fStar != tmpStar1.starNum))
                                            {
                                                tmpStar3.starNum = nCatalog[pm].fStar;
                                                threeClosed = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                                if (tmpStar3.starNum != -1)
                                {
                                    ++cnt;
                                    assingObjectIndexToStar(&candidates[0], &candidates[distIndexFirstStar], sObjects);
                                    tmpStar1.objNum = sObjects[0];
                                    tmpStar2.objNum = sObjects[1];
                                    tmpStar3.objNum = sObjects[2];

                                    for (int c = 0; c < groupCount; c++)
                                    {
                                        if (threeCandidatesCount[c] == 0)
                                        {
                                            threeCandidates[c][0] = tmpStar1;
                                            threeCandidates[c][1] = tmpStar2;
                                            threeCandidates[c][2] = tmpStar3;
                                            threeCandidatesCount[c] = 3;
                                            break;
                                        }
                                        else if (l_st[tmpStar1.starNum] * l_st[threeCandidates[c][0].starNum]
                                                 + m_st[tmpStar1.starNum] * m_st[threeCandidates[c][0].starNum]
                                                 + n_st[tmpStar1.starNum] * n_st[threeCandidates[c][0].starNum] >= CosP)
                                        {
                                            threeCandidates[c][threeCandidatesCount[c]++] = tmpStar1;
                                            threeCandidates[c][threeCandidatesCount[c]++] = tmpStar2;
                                            threeCandidates[c][threeCandidatesCount[c]++] = tmpStar3;
                                            if (threeCandidatesCount[c] > maxGroupSize - 3)
                                            {
                                                qsort(threeCandidates[c], threeCandidatesCount[c], sizeof(StarObject), compareStarObject);
                                                threeCandidatesCount[c] = uniques(threeCandidates[c], threeCandidatesCount[c]);
                                                newPereb(threeCandidates[c], threeCandidatesCount[c], l_st, m_st, n_st);
                                                if (NumDet >= 4)
                                                {
                                                    time2 = rdtsc();
                                                    times.append(time2 - time1);
                                                    //qDebug() /*<< times*/ << NumDet << timer.nsecsElapsed();
                                                    goto next;
                                                }
                                            }


                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // time2 = rdtsc();
        // times.append(time2 - time1);
        //qDebug() << threeCandidates;
        // time1 = rdtsc();
        //qDebug() << groupCount;


        for (int i = 0; i < groupCount; i++)
        {
            if (threeCandidatesCount[i] == 0)
            {
                break;
            }
            qsort(threeCandidates[i], threeCandidatesCount[i], sizeof(StarObject), compareStarObject);
            threeCandidatesCount[i] = uniques(threeCandidates[i], threeCandidatesCount[i]);
            if (threeCandidatesCount[i] >= MinDet)
            {
                //                Mnst = threeCandidatesCount[i] > mPM ? mPM : threeCandidatesCount[i];
                //                for (int i1 = 0; i1 < Mnst - 1; i1++)
                //                {
                //                    for (int j1 = i1 + 1; j1 < Mnst; j1++)
                //                    {
                //                        PM[i1][j1] = PM[j1][i1] = l_st[threeCandidates[i][i1].starNum] * l_st[threeCandidates[i][j1].starNum]
                //                                + m_st[threeCandidates[i][i1].starNum] * m_st[threeCandidates[i][j1].starNum]
                //                                + n_st[threeCandidates[i][i1].starNum] * n_st[threeCandidates[i][j1].starNum];
                //                    }
                //                }
                NumDet = 0;
                //Pereb();
                newPereb(threeCandidates[i], threeCandidatesCount[i], l_st, m_st, n_st);
                if (NumDet >= MinDet)
                {
                    count += NumDet;
                    break;
                }
                NumDet = 0;
            }
        }
        time2 = rdtsc();
        times.append(time2 - time1);
                next:
        qDebug() /*<< times*/ << NumDet << timer.nsecsElapsed();

    }


    qDebug() << "FFFF" << count;
}


void MainWindow::IdentObjects(unsigned short graf[MaxDef][MaxDef], unsigned short Nf, unsigned short* im_save, unsigned short* Nmax)
{
    int i, j, N, min, i_del, num_del, k, buf, fl;
    float eps_del, eps_cur;    //16/11/2016//

    for (i=0;i<Nf;i++)
    {
        im[i]=i;
    }

    *Nmax=0;
    num_del=0;
    N=Nf;

    while (*Nmax<N)
    {
        num_del=0;
        while(N>1) //??? MinDet - ?
        {
            //iteration
            min=1000;
            for (i=0;i<N;i++)
            {
                sum_row[im[i]]=0;
                for (j=0;j<N;j++)
                    sum_row[im[i]]+=graf[im[i]][im[j]];
                if (sum_row[im[i]]<=min)
                {
                    min =sum_row[im[i]];
                    i_del=im[i];
                }
            }
            if (N==min) break;
            ///////////////16/11/2016//////////////////////////////
            for (i=0;i<N;i++)
            {
                if ((sum_row[im[i]]==min)&&(im[i]!=i_del))
                {
                    eps_del=0; eps_cur=0;
                    for (k=0; k<N; k++)
                    {
                        eps_del+=MatDif[i_del][im[k]];
                        eps_cur+=MatDif[im[i]][im[k]];
                    }
                    if (eps_cur>eps_del)
                    {
                        min =sum_row[im[i]];
                        i_del=im[i];
                    }
                }
            }
            ///////////////////////////////////////////////////////


            for(i=0;i<N;i++)
                if (im[i]>i_del)
                    im[i-1]=im[i];

            im_del[num_del]=i_del;
            num_del++;
            N--;
        }  //end while (N>1)

        if (N>*Nmax)
        {
            for(i=0;i<N;i++)
                im_save[i]=im[i];
            *Nmax=N;
        }

        if (num_del>*Nmax)  //16/11/2016//
        {
            fl=1;
            while(fl)
            {
                fl=0;
                k=0;
                while (k<num_del-1)
                {
                    if (im_del[k+1]<im_del[k])
                    {
                        buf=im_del[k];
                        im_del[k]=im_del[k+1];
                        im_del[k+1]=buf; fl=1;
                    }
                    k++;
                }
            }

            for (i=0;i<num_del;i++)
                im[i]=im_del[i];
        }
        N=num_del;
    }  //end while (Nmax<N)
}


void MainWindow::newPereb(StarObject* threeCandidates, int count, QVector <double>& l_st, QVector <double>& m_st, QVector <double>& n_st)
{
    for (int i = 0; i < count; i++)
    {
        for (int j = 0; j < count; j++)
        {
            MatCmp[i][j] = 1;
        }
    }
    float Lf, Mf, Nf, Ls, Ms, Ns;
    for (int i = 0; i < count - 1; i++)
    {
        Lf = l_st[threeCandidates[i].starNum];
        Mf = m_st[threeCandidates[i].starNum];
        Nf = n_st[threeCandidates[i].starNum];

        for (int j = i + 1; j < count; j++)
        {
            Ls = l_st[threeCandidates[j].starNum];
            Ms = m_st[threeCandidates[j].starNum];
            Ns = n_st[threeCandidates[j].starNum];
            if((MatDif[i][j] = fabs((Ncos[threeCandidates[i].objNum][threeCandidates[j].objNum]) - (Lf * Ls + Mf * Ms + Nf * Ns)))
                    > Eps[threeCandidates[i].objNum][threeCandidates[j].objNum])
            {
                MatCmp[i][j] = 0;
                MatCmp[j][i] = 0;
            }
            MatDif[j][i] = MatDif[i][j];
        }
    }
    IdentObjects(MatCmp, count, NumIdent, &NumDet);
}




