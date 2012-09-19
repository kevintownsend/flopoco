#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

int main () {


  string target="Virtex5";
  int wX=1;
  int wY=1;
  int wOut=2;
  bool signedIO=0;
  float ratio=1.0;
  bool enableSuperTiles=1;
 int smth=1;
 int smth1=1;

  ofstream myfile;
  printf("This will take a time... be patient \n");
  /*******************************************************************************************************************************************/
  myfile.open ("Virtex5.cmd");
  myfile << "#IntMultiplier\n";
  target="Virtex5";

  printf("Generating file for Virtex5 \n");
  for(wX=1;wX<=60;wX=wX+3)
    for(wY=1;wY<=60;wY=wY+5)
        for(wOut=min(wX,wY);wOut<=wX+wY;wOut=wOut+10)
          {
            if(smth1==1)
            {  smth1=0;
               signedIO=1;
            }
            else
            {   smth1=1;
                signedIO=0;
            }

             for(ratio=0.0;ratio<=1.0;ratio=ratio+0.5)

                {
                    if(smth==1)
                    {
                        enableSuperTiles=1;
                        smth=0;
                    }
                    else
                    {
                        enableSuperTiles=0;
                        smth=1;
                    }
                    myfile<<"flopoco -target="<<target<<" IntMultiplier "<<wX<<" "<<wY<<" "<<wOut<<" "<<signedIO<<" "<<ratio<<" "<<enableSuperTiles<<"\n";
                }
          }

  myfile.close();

/*******************************************************************************************************************************************/


/*******************************************************************************************************************************************/
  myfile.open ("Virtex4.cmd");
  myfile << "#IntMultiplier\n";
  target="Virtex4";
   printf("Generating file for Virtex4 \n");
   for(wX=1;wX<=60;wX=wX+3)
    for(wY=1;wY<=60;wY=wY+5)
        for(wOut=min(wX,wY);wOut<=wX+wY;wOut=wOut+10)
          {
            if(smth1==1)
            {  smth1=0;
               signedIO=1;
            }
            else
            {   smth1=1;
                signedIO=0;
            }

             for(ratio=0.0;ratio<=1.0;ratio=ratio+0.5)

                {
                    if(smth==1)
                    {
                        enableSuperTiles=1;
                        smth=0;
                    }
                    else
                    {
                        enableSuperTiles=0;
                        smth=1;
                    }
                    myfile<<"flopoco -target="<<target<<" IntMultiplier "<<wX<<" "<<wY<<" "<<wOut<<" "<<signedIO<<" "<<ratio<<" "<<enableSuperTiles<<"\n";
                }
          }
  myfile.close();
/*******************************************************************************************************************************************/


/*******************************************************************************************************************************************/
  myfile.open ("Spartan3.cmd");
  myfile << "#IntMultiplier\n";
  target="Spartan3";
   printf("Generating file for Spartan3 \n");
   for(wX=1;wX<=60;wX=wX+3)
    for(wY=1;wY<=60;wY=wY+5)
        for(wOut=min(wX,wY);wOut<=wX+wY;wOut=wOut+10)
          {
            if(smth1==1)
            {  smth1=0;
               signedIO=1;
            }
            else
            {   smth1=1;
                signedIO=0;
            }

             for(ratio=0.0;ratio<=1.0;ratio=ratio+0.5)

                {
                    if(smth==1)
                    {
                        enableSuperTiles=1;
                        smth=0;
                    }
                    else
                    {
                        enableSuperTiles=0;
                        smth=1;
                    }
                    myfile<<"flopoco -target="<<target<<" IntMultiplier "<<wX<<" "<<wY<<" "<<wOut<<" "<<signedIO<<" "<<ratio<<" "<<enableSuperTiles<<"\n";
                }
          }
  myfile.close();
  /*******************************************************************************************************************************************/
  printf("\n\n\n*********************************************************************");
  printf("Xilinx files have been generated \n");
  printf("Altera files are generating now \n");
  printf("\n\n\n*********************************************************************");

/*******************************************************************************************************************************************/
  myfile.open ("StratixII.cmd");
  myfile << "#IntMultiplier\n";
  target="StratixII";
  printf("Generating file for Stratix2 \n");
 for(wX=1;wX<=60;wX=wX+3)
    for(wY=1;wY<=60;wY=wY+5)
        for(wOut=min(wX,wY);wOut<=wX+wY;wOut=wOut+10)
          {
            if(smth1==1)
            {  smth1=0;
               signedIO=1;
            }
            else
            {   smth1=1;
                signedIO=0;
            }

             for(ratio=0.0;ratio<=1.0;ratio=ratio+0.5)

                {
                    if(smth==1)
                    {
                        enableSuperTiles=1;
                        smth=0;
                    }
                    else
                    {
                        enableSuperTiles=0;
                        smth=1;
                    }
                    myfile<<"flopoco -target="<<target<<" IntMultiplier "<<wX<<" "<<wY<<" "<<wOut<<" "<<signedIO<<" "<<ratio<<" "<<enableSuperTiles<<"\n";
                }
          }

  myfile.close();
  /*******************************************************************************************************************************************/


/*******************************************************************************************************************************************/
  myfile.open ("StratixIII.cmd");
  myfile << "#IntMultiplier\n";
   printf("Generating file for Stratix3 \n");
  target="StratixIII";
   for(wX=1;wX<=60;wX=wX+3)
    for(wY=1;wY<=60;wY=wY+5)
        for(wOut=min(wX,wY);wOut<=wX+wY;wOut=wOut+10)
          {
            if(smth1==1)
            {  smth1=0;
               signedIO=1;
            }
            else
            {   smth1=1;
                signedIO=0;
            }

             for(ratio=0.0;ratio<=1.0;ratio=ratio+0.5)

                {
                    if(smth==1)
                    {
                        enableSuperTiles=1;
                        smth=0;
                    }
                    else
                    {
                        enableSuperTiles=0;
                        smth=1;
                    }
                    myfile<<"flopoco -target="<<target<<" IntMultiplier "<<wX<<" "<<wY<<" "<<wOut<<" "<<signedIO<<" "<<ratio<<" "<<enableSuperTiles<<"\n";
                }
          }
          
     myfile.close();     
     
     
  myfile.open ("StratixIV.cmd");
  myfile << "#IntMultiplier\n";
   printf("Generating file for Stratix4 \n");
  target="StratixIV";
   for(wX=1;wX<=60;wX=wX+3)
    for(wY=1;wY<=60;wY=wY+5)
        for(wOut=min(wX,wY);wOut<=wX+wY;wOut=wOut+10)
          {
            if(smth1==1)
            {  smth1=0;
               signedIO=1;
            }
            else
            {   smth1=1;
                signedIO=0;
            }

             for(ratio=0.0;ratio<=1.0;ratio=ratio+0.5)

                {
                    if(smth==1)
                    {
                        enableSuperTiles=1;
                        smth=0;
                    }
                    else
                    {
                        enableSuperTiles=0;
                        smth=1;
                    }
                    myfile<<"flopoco -target="<<target<<" IntMultiplier "<<wX<<" "<<wY<<" "<<wOut<<" "<<signedIO<<" "<<ratio<<" "<<enableSuperTiles<<"\n";
                }
          }
          
     myfile.close();     

printf("\n\n\n*********************************************************************");
printf("\n\n\nGENERATION COMPLETED");
printf("\n\n\n*********************************************************************");
return 1;
}

