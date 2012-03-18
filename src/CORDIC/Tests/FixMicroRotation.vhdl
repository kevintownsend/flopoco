--------------------------------------------------------------------------------
--                            IntAdder_6_f400_uid4
--                      (IntAdderClassical_6_f400_uid6)
-- This operator is part of the Infinite Virtual Library FloPoCoLib
-- All rights reserved 
-- Authors: Bogdan Pasca, Florent de Dinechin (2008-2010)
--------------------------------------------------------------------------------
-- combinatorial

library ieee;
use ieee.std_logic_1164.all;
use ieee.std_logic_arith.all;
use ieee.std_logic_unsigned.all;
library std;
use std.textio.all;
library work;

entity IntAdder_6_f400_uid4 is
   port ( X : in  std_logic_vector(5 downto 0);
          Y : in  std_logic_vector(5 downto 0);
          Cin : in std_logic;
          R : out  std_logic_vector(5 downto 0)   );
end entity;

architecture arch of IntAdder_6_f400_uid4 is
begin
   --Classical
    R <= X + Y + Cin;
end architecture;

--------------------------------------------------------------------------------
--                       FixMicroRotation_6_uid2_stage0
-- This operator is part of the Infinite Virtual Library FloPoCoLib
-- All rights reserved 
-- Authors: Istoan Matei, Florent de Dinechin (2008-2012)
--------------------------------------------------------------------------------
-- combinatorial

library ieee;
use ieee.std_logic_1164.all;
use ieee.std_logic_arith.all;
use ieee.std_logic_unsigned.all;
library std;
use std.textio.all;
library work;

entity FixMicroRotation_6_uid2_stage0 is
   port ( Xin : in  std_logic_vector(5 downto 0);
          Yin : in  std_logic_vector(5 downto 0);
          Zin : in  std_logic_vector(5 downto 0);
          Din : in std_logic;
          Xout : out  std_logic_vector(5 downto 0);
          Yout : out  std_logic_vector(5 downto 0);
          Zout : out  std_logic_vector(5 downto 0);
          Dout : out std_logic   );
end entity;

architecture arch of FixMicroRotation_6_uid2_stage0 is
   component IntAdder_6_f400_uid4 is
      port ( X : in  std_logic_vector(5 downto 0);
             Y : in  std_logic_vector(5 downto 0);
             Cin : in std_logic;
             R : out  std_logic_vector(5 downto 0)   );
   end component;

signal XinShift :  std_logic_vector(5 downto 0);
signal YinShift :  std_logic_vector(5 downto 0);
signal newXinShift :  std_logic_vector(5 downto 0);
signal newYinShift :  std_logic_vector(5 downto 0);
signal cInNewX : std_logic;
signal cInNewY : std_logic;
signal intXout :  std_logic_vector(5 downto 0);
signal intYout :  std_logic_vector(5 downto 0);
signal atan2PowStage :  std_logic_vector(5 downto 0);
signal newAtan2PowStage :  std_logic_vector(5 downto 0);
signal cInZ : std_logic;
signal intZout :  std_logic_vector(5 downto 0);
signal intDout : std_logic;
begin
   XinShift <= Xin;
   YinShift <= Yin;
   newXinShift <= XinShift xor (5 downto 0 => Din);
   newYinShift <= YinShift xor (5 downto 0 => (not Din));
   cInNewX<= Din;
   cInNewY<= not Din;
   xAdder: IntAdder_6_f400_uid4
      port map ( Cin => cInNewY,
                 R => intXout,
                 X => Xin,
                 Y => newYinShift);

   yAdder: IntAdder_6_f400_uid4
      port map ( Cin => cInNewX,
                 R => intYout,
                 X => Yin,
                 Y => newXinShift);

   atan2PowStage <= '0' & "00011";
   newAtan2PowStage <= atan2PowStage xor (5 downto 0 => (not Din));
   cInZ<= not Din;
   zAdder: IntAdder_6_f400_uid4
      port map ( Cin => cInZ,
                 R => intZout,
                 X => Zin,
                 Y => newAtan2PowStage);

   intDout <= intZout(5);
   Xout <= intXout;
   Yout <= intYout;
   Zout <= intZout;
   Dout <= intDout;
end architecture;

--------------------------------------------------------------------------------
--                  TestBench_FixMicroRotation_6_uid2_stage0
-- This operator is part of the Infinite Virtual Library FloPoCoLib
-- All rights reserved 
-- Authors: Florent de Dinechin, Cristian Klein (2007)
--------------------------------------------------------------------------------
library ieee;
use ieee.std_logic_1164.all;
use ieee.std_logic_arith.all;
use ieee.std_logic_unsigned.all;
library std;
use std.textio.all;
library work;

entity TestBench_FixMicroRotation_6_uid2_stage0 is
end entity;

architecture behavorial of TestBench_FixMicroRotation_6_uid2_stage0 is
   component FixMicroRotation_6_uid2_stage0 is
      port ( Xin : in  std_logic_vector(5 downto 0);
             Yin : in  std_logic_vector(5 downto 0);
             Zin : in  std_logic_vector(5 downto 0);
             Din : in std_logic;
             Xout : out  std_logic_vector(5 downto 0);
             Yout : out  std_logic_vector(5 downto 0);
             Zout : out  std_logic_vector(5 downto 0);
             Dout : out std_logic   );
   end component;
   signal Xin :  std_logic_vector(5 downto 0);
   signal Yin :  std_logic_vector(5 downto 0);
   signal Zin :  std_logic_vector(5 downto 0);
   signal Din : std_logic;
   signal Xout :  std_logic_vector(5 downto 0);
   signal Yout :  std_logic_vector(5 downto 0);
   signal Zout :  std_logic_vector(5 downto 0);
   signal Dout : std_logic;
   signal clk : std_logic;
   signal rst : std_logic;

   -- FP compare function (found vs. real)
   function fp_equal(a : std_logic_vector; b : std_logic_vector) return boolean is
   begin
      if b(b'high downto b'high-1) = "01" then
         return a = b;
      elsif b(b'high downto b'high-1) = "11" then
         return (a(a'high downto a'high-1)=b(b'high downto b'high-1));
      else
         return a(a'high downto a'high-2) = b(b'high downto b'high-2);
      end if;
   end;



 -- converts std_logic into a character
   function chr(sl: std_logic) return character is
      variable c: character;
   begin
      case sl is
         when 'U' => c:= 'U';
         when 'X' => c:= 'X';
         when '0' => c:= '0';
         when '1' => c:= '1';
         when 'Z' => c:= 'Z';
         when 'W' => c:= 'W';
         when 'L' => c:= 'L';
         when 'H' => c:= 'H';
         when '-' => c:= '-';
      end case;
      return c;
   end chr;
   -- converts bit to std_logic (1 to 1)
   function to_stdlogic(b : bit) return std_logic is
       variable sl : std_logic;
   begin
      case b is 
         when '0' => sl := '0';
         when '1' => sl := '1';
      end case;
      return sl;
   end to_stdlogic;
   -- converts std_logic into a string (1 to 1)
   function str(sl: std_logic) return string is
    variable s: string(1 to 1);
    begin
      s(1) := chr(sl);
      return s;
   end str;
   -- converts std_logic_vector into a string (binary base)
   -- (this also takes care of the fact that the range of
   --  a string is natural while a std_logic_vector may
   --  have an integer range)
   function str(slv: std_logic_vector) return string is
      variable result : string (1 to slv'length);
      variable r : integer;
   begin
      r := 1;
      for i in slv'range loop
         result(r) := chr(slv(i));
         r := r + 1;
      end loop;
      return result;
   end str;




   -- test isZero
   function iszero(a : std_logic_vector) return boolean is
   begin
      return  a = (a'high downto 0 => '0');
   end;


   -- FP IEEE compare function (found vs. real)
   function fp_equal_ieee(a : std_logic_vector; b : std_logic_vector; we : integer; wf : integer) return boolean is
   begin
      if a(wf+we downto wf) = b(wf+we downto wf) and b(we+wf-1 downto wf) = (we downto 1 => '1') then
         if iszero(b(wf-1 downto 0)) then return  iszero(a(wf-1 downto 0));
         else return not iszero(a(wf - 1 downto 0));
         end if;
      else
         return a(a'high downto 0) = b(b'high downto 0);
      end if;
   end;
begin
   test: FixMicroRotation_6_uid2_stage0
      port map ( Din => Din,
                 Dout => Dout,
                 Xin => Xin,
                 Xout => Xout,
                 Yin => Yin,
                 Yout => Yout,
                 Zin => Zin,
                 Zout => Zout);
   -- Ticking clock signal
   process
   begin
      clk <= '0';
      wait for 5 ns;
      clk <= '1';
      wait for 5 ns;
   end process;

   -- Reading the input from a file 
   process
      variable inline : line; 
      variable counter : integer := 1;
      variable errorCounter : integer := 0;
      variable possibilityNumber : integer := 0;
      variable localErrorCounter : integer := 0;
      variable tmpChar : character;
      file inputsFile : text is "test.input"; 
      variable V_Xin : bit_vector(5 downto 0);
      variable V_Yin : bit_vector(5 downto 0);
      variable V_Zin : bit_vector(5 downto 0);
      variable V_Din : bit_vector(0 downto 0);
      variable V_Xout : bit_vector(5 downto 0);
      variable V_Yout : bit_vector(5 downto 0);
      variable V_Zout : bit_vector(5 downto 0);
      variable V_Dout : bit_vector(0 downto 0);
   begin
      -- Send reset
      rst <= '1';
      wait for 10 ns;
      rst <= '0';
      while not endfile(inputsFile) loop
          -- positionning inputs
         readline(inputsFile,inline);
         read(inline ,V_Xin);
         read(inline,tmpChar);
         Xin <= to_stdlogicvector(V_Xin);
         read(inline ,V_Yin);
         read(inline,tmpChar);
         Yin <= to_stdlogicvector(V_Yin);
         read(inline ,V_Zin);
         read(inline,tmpChar);
         Zin <= to_stdlogicvector(V_Zin);
         read(inline ,V_Din);
         read(inline,tmpChar);
         Din <= to_stdlogicvector(V_Din)(0);
         readline(inputsFile,inline);
         wait for 10 ns;
      end loop;
      wait for 10000 ns; -- wait for simulation to finish
   end process;
          -- verifying the corresponding output
         process
      variable inline : line; 
      variable inline0 : line; 
      variable counter : integer := 1;
      variable errorCounter : integer := 0;
      variable possibilityNumber : integer := 0;
      variable localErrorCounter : integer := 0;
      variable tmpChar : character;
      file inputsFile : text is "test.input"; 
      variable V_Xin : bit_vector(5 downto 0);
      variable V_Yin : bit_vector(5 downto 0);
      variable V_Zin : bit_vector(5 downto 0);
      variable V_Din : bit;
      variable V_Xout : bit_vector(5 downto 0);
      variable V_Yout : bit_vector(5 downto 0);
      variable V_Zout : bit_vector(5 downto 0);
      variable V_Dout : bit;
   begin
          wait for 10 ns;
      wait for 2 ns; -- wait for pipeline to flush
      while not endfile(inputsFile) loop
          -- positionning inputs
         readline(inputsFile,inline0);
         readline(inputsFile,inline);
         read(inline, possibilityNumber);
         localErrorCounter := 0;
         read(inline,tmpChar);
         if possibilityNumber = 0 then
            localErrorCounter := 0;
         elsif possibilityNumber = 1 then 
            read(inline ,V_Xout);
            if not (Xout= to_stdlogicvector(V_Xout)) then 
               assert false report("Incorrect output for Xout,expected value: " & str(to_stdlogicvector(V_Xout)) & ", result: " & str(Xout)) &  "|| line : " & integer'image(counter) & " of input file " ;
                errorCounter := errorCounter + 1;
            end if;
         else
            for i in possibilityNumber downto 1 loop 
               read(inline ,V_Xout);
               read(inline,tmpChar);
               if (Xout= to_stdlogicvector(V_Xout))  then localErrorCounter := 1; end if;
            end loop;
             if (localErrorCounter = 0) then 
               errorCounter := errorCounter + 1; -- incrementing global error counter
               assert false report("Incorrect output for Xout, expected value : " & str(to_stdlogicvector(V_Xout)) & "... (other values line " & integer'image(counter) & " of test.input), result:  " & str(Xout) &  "|| line : " & integer'image(counter) & " of input file ") ;
            end if;
         end if;
         read(inline, possibilityNumber);
         localErrorCounter := 0;
         read(inline,tmpChar);
         if possibilityNumber = 0 then
            localErrorCounter := 0;
         elsif possibilityNumber = 1 then 
            read(inline ,V_Yout);
            if not (Yout= to_stdlogicvector(V_Yout)) then 
               assert false report("Incorrect output for Yout,expected value: " & str(to_stdlogicvector(V_Yout)) & ", result: " & str(Yout)) &  "|| line : " & integer'image(counter) & " of input file " ;
                errorCounter := errorCounter + 1;
            end if;
         else
            for i in possibilityNumber downto 1 loop 
               read(inline ,V_Yout);
               read(inline,tmpChar);
               if (Yout= to_stdlogicvector(V_Yout))  then localErrorCounter := 1; end if;
            end loop;
             if (localErrorCounter = 0) then 
               errorCounter := errorCounter + 1; -- incrementing global error counter
               assert false report("Incorrect output for Yout, expected value : " & str(to_stdlogicvector(V_Yout)) & "... (other values line " & integer'image(counter) & " of test.input), result:  " & str(Yout) &  "|| line : " & integer'image(counter) & " of input file ") ;
            end if;
         end if;
         read(inline, possibilityNumber);
         localErrorCounter := 0;
         read(inline,tmpChar);
         if possibilityNumber = 0 then
            localErrorCounter := 0;
         elsif possibilityNumber = 1 then 
            read(inline ,V_Zout);
            if not (Zout= to_stdlogicvector(V_Zout)) then 
               assert false report("Incorrect output for Zout,expected value: " & str(to_stdlogicvector(V_Zout)) & ", result: " & str(Zout)) &  "|| line : " & integer'image(counter) & " of input file " ;
                errorCounter := errorCounter + 1;
            end if;
         else
            for i in possibilityNumber downto 1 loop 
               read(inline ,V_Zout);
               read(inline,tmpChar);
               if (Zout= to_stdlogicvector(V_Zout))  then localErrorCounter := 1; end if;
            end loop;
             if (localErrorCounter = 0) then 
               errorCounter := errorCounter + 1; -- incrementing global error counter
               assert false report("Incorrect output for Zout, expected value : " & str(to_stdlogicvector(V_Zout)) & "... (other values line " & integer'image(counter) & " of test.input), result:  " & str(Zout) &  "|| line : " & integer'image(counter) & " of input file ") ;
            end if;
         end if;
         read(inline, possibilityNumber);
         localErrorCounter := 0;
         read(inline,tmpChar);
         if possibilityNumber = 0 then
            localErrorCounter := 0;
         elsif possibilityNumber = 1 then 
            read(inline ,V_Dout);
            if not (Dout= to_stdlogic(V_Dout)) then 
               assert false report("Incorrect output for Dout,expected value: " & str(to_stdlogic(V_Dout)) & ", result: " & str(Dout)) &  "|| line : " & integer'image(counter) & " of input file " ;
                errorCounter := errorCounter + 1;
            end if;
         else
            for i in possibilityNumber downto 1 loop 
               read(inline ,V_Dout);
               read(inline,tmpChar);
               if (Dout= to_stdlogic(V_Dout))  then localErrorCounter := 1; end if;
            end loop;
             if (localErrorCounter = 0) then 
               errorCounter := errorCounter + 1; -- incrementing global error counter
               assert false report("Incorrect output for Dout, expected value: " & str(to_stdlogic(V_Dout)) & "... (other values line " & integer'image(counter) & " of test.input), result:  " & str(Dout) &  "|| line : " & integer'image(counter) & " of input file ") ;
            end if;
         end if;
          wait for 10 ns; -- wait for pipeline to flush
         counter := counter + 2;
      end loop;
      report (integer'image(errorCounter) & " error(s) encoutered.");
      assert false report "End of simulation" severity failure;
   end process;

end architecture;

