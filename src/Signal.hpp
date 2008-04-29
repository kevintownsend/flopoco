#ifndef __SIGNAL_HPP
#define __SIGNAL_HPP

#include <iostream>
#include <sstream>

// TODO: Cleanup & Doxygenize.

	// A local class to be replaced by the cpphdl version some day
	class Signal{
	public:
		typedef enum{in,out,wire,registered,registered_with_async_reset,registered_with_sync_reset } type_t;
		
		Signal(const std::string name, const type_t type, const int width = 1): 
			_name(name), _type(type), _wE(0), _wF(0), _width(width){
			_isFP = false;
		}

		Signal(const std::string name, const type_t type, const int wE, const int wF): 
			_name(name), _type(type), _wE(wE), _wF(wF), _width(wE+wF+3){
			_isFP = true;
		}
		
		~Signal(){}
		
		std::string id(){return _name;}
		
		int width(){return _width;}
		int wE(){return(_wE);}
		int wF(){return(_wF);}
		bool isFP(){return _isFP;}
		
		type_t type() {return _type;}
		
		std::string toVHDL() {
			std::ostringstream o; 
			if(type()==Signal::wire || type()==Signal::registered || type()==Signal::registered_with_async_reset || type()==Signal::registered_with_sync_reset) 
				o << "signal ";
			o << id();
//       if(type()==Signal::registered || type()==Signal::registered_with_reset|| type()==Signal::delay)
// 	o << ", " << id() << "_d";	 
			o << " : ";
			if(type()==Signal::in)
				 o << "in ";
			if(type()==Signal::out)
				 o << "out ";

			if (1==width()) 
				o << "std_logic" ;
			else 
				if(_isFP) 
					o << " std_logic_vector(" << wE() <<"+"<<wF() << "+2 downto 0)";
				else
					o << " std_logic_vector(" << width()-1 << " downto 0)";
			return o.str();
		}

	private:
		const std::string _name;
		const type_t _type;
		const uint32_t _width;
		bool _isFP;
		const uint32_t _wE;  // used only for FP signals
		const uint32_t _wF;  // used only for FP signals
	};

#endif

