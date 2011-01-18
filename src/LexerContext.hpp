#ifndef __LEXER_CONTEXT_HPP__
#define __LEXER_CONTEXT_HPP__

#include <iostream>

using namespace std;

class LexerContext {
	public:
		void* scanner;
		int result;
		istream* is;
		ostream* os;

	public:
		LexerContext(istream* is = &cin, ostream* os = &cout) {
			init_scanner();
			this->is = is;
			this->os = os;
		}

		void lex();

		virtual ~LexerContext() { destroy_scanner();}

	protected:
		void init_scanner();

		void destroy_scanner();
};



#endif
