#include "fluxes.h"
#include "visc.h"
#include "ideal.h"
#include "isotherm.h"
#include "loadData.h"

int main()
{
    Mat QOne = loadFromTxt("startingDistributions/2DShocktubeQ1.txt");
    Mat QTwox = loadFromTxt("startingDistributions/2DShocktubeQ2x.txt");
    Mat QTwoy = loadFromTxt("startingDistributions/2DShocktubeQ2y.txt");
    Mat QThree = loadFromTxt("startingDistributions/2DShocktubeQ3.txt");
    viscSimulation A(QOne, QTwox, QTwoy, QThree, 0.01, 0.01, 3.0);
    //A.print();
    A.update(50);

    A.saveTo("data.txt");
}
