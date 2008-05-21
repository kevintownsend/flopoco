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
			_name(name), _type(type), _wE(0), _wF(0), _width(width), _low(0), _high(width-1), _isSubSignal(false) {
			_isFP = false;
		}

		Signal(const std::string name, const type_t type, const int wE, const int wF): 
			_name(name), _type(type), _wE(wE), _wF(wF), _width(wE+wF+3), _low(0), _high(_width-1), _isSubSignal(false) {
			_isFP = true;
		}
		
		~Signal(){}
		
		std::string id() {
			if (_isSubSignal == false)
				return _name;
			else
			{
				std::stringstream o;
				if (_width == 1)
					o << _name << "(" << _low << ")";
				else
					o << _name << "(" << _high << " downto " << _low << ")";
				return o.str();
			}
		}
		
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

		Signal getSubSignal(int low, int high)
		{
			if (low < _low)
				throw std::string("Attempted to return subsignal with smaller low index.");
			if (high > _high)
				throw std::string("Attempted to return subsignal with bigger high index."); 

			Signal s(_name, _type, high-low+1);
			s._low = low;
			s._high = high;
			s._isSubSignal = true;

			return s;
		}

		Signal getException()
		{
			if (!_isFP)
				throw std::string("Not a floating point signal.");
			return getSubSignal(1+_wE+_wF, 2+_wE+_wF);
		}

		Signal getSign()
		{
			if (!_isFP)
				throw std::string("Not a floating point signal.");
			return getSubSignal(_wE+_wF, _wE+_wF);
		}

		Signal getExponent()
		{
			if (!_isFP)
				throw std::string("Not a floating point signal.");
			return getSubSignal(_wF, _wE+_wF-1);
		}

		Signal getMantissa()
		{
			if (!_isFP)
				throw std::string("Not a floating point signal.");
			return getSubSignal(0, _wF-1);
		}

	private:
		std::string _name;
		type_t _type;
		uint32_t _width;
		bool _isFP;
		uint32_t _wE;  // used only for FP signals
		uint32_t _wF;  // used only for FP signals
		bool _isSubSignal;
		uint32_t _low, _high; // used only for „subsignals”, which are non FP
	};

#endif

