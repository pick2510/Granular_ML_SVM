#ifndef INTERFACES_INCLUDED
#define INTERFACES_INCLUDED
#include <iostream>
#include <string>
#include "PeriodicCellList.h"
#include "Plots.h"
int InvStatMechCLI();
int PairStatisticsCLI();
int CollectiveCoordinateCLI();

//defined in debug.cpp
int Debug();
//defined in debug.cpp
int CalculateTau(void);
int CalculateEvAndQuantizerError();
void FixDe();

int cli();
#endif