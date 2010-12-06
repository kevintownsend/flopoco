%{
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include "ExpressionParser.h"

int FlopocoExpressionlex(void);
int FlopocoExpressionparse(void);
int FlopocoExpressionerror(char *);

using namespace std;
%}

A [/]
B [*]
C [^/*]

%option noyywrap
%%

{A}({A}.*|{B}({A}|{C}|{B}+{C})*{B}+{A}) { /* skip comment */ }

" "               {
                  }
"\t"              {
                  }
                  
                  
"+"               {
                    	FlopocoExpressionlval.c_type = FlopocoExpressiontext[0];
						cout << " Found +" << endl;
	                    return FPEXPRESSION_PLUS;
                  }
"-"               {
					cout << " Found -" << endl;
					FlopocoExpressionlval.c_type = FlopocoExpressiontext[0];
                   return FPEXPRESSION_MINUS;
                  }
"*"               { 
					cout << " Found *" << endl;
					return FPEXPRESSION_TIMES;}
"/"               { 
					cout << " Found /" << endl;
					return FPEXPRESSION_DIV;}                  
                  
"="               {
					cout << "Found assignment" << endl;
                    FlopocoExpressionlval.c_type = FlopocoExpressiontext[0];
                    return FPEXPRESSION_EQUALS;
                  }

";"               { 
					return FPEXPRESSION_PV;
				  }

"output"          { 
					cout << "Found output keyword" << endl;
					return FPEXPRESSION_OUTPUT;
				  }

"sqrt"            {
					cout << "Found sqrt" << endl;
					return FPEXPRESSION_SQRT;
					}             
"sqr"            {
					cout << "Found sqr" << endl;
					return FPEXPRESSION_SQR;
				}
"exp"            {
					cout << "Found exp" << endl;
					return FPEXPRESSION_EXP;
					}             
"log"            {
					cout << "Found log" << endl;
					return FPEXPRESSION_LOG;
				}
								  
"("				  { cout << "Found (" << endl; return FPEXPRESSION_LPAR;}
")"				  { cout << "Found )" << endl; return FPEXPRESSION_RPAR;}


[a-zA-Z]+            {
						
						FlopocoExpressionlval.s_type = (char*) malloc( sizeof(char)*strlen(FlopocoExpressiontext) );
						for (unsigned i=0;i<strlen(FlopocoExpressiontext);i++)
							FlopocoExpressionlval.s_type[i] = FlopocoExpressiontext[i];
						cout << "found variable with name: " << FlopocoExpressionlval.s_type << endl;
						return FPEXPRESSION_VARIABLE;
				  }

[1-9][[0-9]*  {
                    FlopocoExpressionlval.i_type = atoi(FlopocoExpressiontext);
                    cout << "found number " << FlopocoExpressionlval.i_type << endl;
                    return FPEXPRESSION_NUMBER;
                  }
.                 {
                    FlopocoExpressionlval.c_type = FlopocoExpressiontext[0];
                    return FPEXPRESSION_UNKNOWN;
                  }
<<EOF>>           {
                    yyterminate();
                  }

%%

/*
int main(void) {
    yylex();
    return 0;
} */

