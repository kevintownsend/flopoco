--------------------------------------------------------------------------------
--                            IntAdder_8_f400_uid6
--                      (IntAdderClassical_8_f400_uid8)
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

entity IntAdder_8_f400_uid6 is
   port ( X : in  std_logic_vector(7 downto 0);
          Y : in  std_logic_vector(7 downto 0);
          Cin : in std_logic;
          R : out  std_logic_vector(7 downto 0)   );
end entity;

architecture arch of IntAdder_8_f400_uid6 is
begin
   --Classical
    R <= X + Y + Cin;
end architecture;

--------------------------------------------------------------------------------
--                       FixMicroRotation_8_uid4_stage0
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

entity FixMicroRotation_8_uid4_stage0 is
   port ( Xin : in  std_logic_vector(7 downto 0);
          Yin : in  std_logic_vector(7 downto 0);
          Zin : in  std_logic_vector(7 downto 0);
          Din : in std_logic;
          Xout : out  std_logic_vector(7 downto 0);
          Yout : out  std_logic_vector(7 downto 0);
          Zout : out  std_logic_vector(7 downto 0);
          Dout : out std_logic   );
end entity;

architecture arch of FixMicroRotation_8_uid4_stage0 is
   component IntAdder_8_f400_uid6 is
      port ( X : in  std_logic_vector(7 downto 0);
             Y : in  std_logic_vector(7 downto 0);
             Cin : in std_logic;
             R : out  std_logic_vector(7 downto 0)   );
   end component;

signal XinShift :  std_logic_vector(7 downto 0);
signal YinShift :  std_logic_vector(7 downto 0);
signal newXinShift :  std_logic_vector(7 downto 0);
signal newYinShift :  std_logic_vector(7 downto 0);
signal cInNewX : std_logic;
signal cInNewY : std_logic;
signal intXout :  std_logic_vector(7 downto 0);
signal intYout :  std_logic_vector(7 downto 0);
signal atan2PowStage :  std_logic_vector(7 downto 0);
signal newAtan2PowStage :  std_logic_vector(7 downto 0);
signal cInZ : std_logic;
signal intZout :  std_logic_vector(7 downto 0);
signal intDout : std_logic;
begin
   XinShift <= Xin;
   YinShift <= Yin;
   newXinShift <= XinShift xor (7 downto 0 => Din);
   newYinShift <= YinShift xor (7 downto 0 => (not Din));
   cInNewX<= Din;
   cInNewY<= not Din;
   xAdder: IntAdder_8_f400_uid6
      port map ( Cin => cInNewY,
                 R => intXout,
                 X => Xin,
                 Y => newYinShift);

   yAdder: IntAdder_8_f400_uid6
      port map ( Cin => cInNewX,
                 R => intYout,
                 X => Yin,
                 Y => newXinShift);

   atan2PowStage <= '0' & "0001101";
   newAtan2PowStage <= atan2PowStage xor (7 downto 0 => (not Din));
   cInZ<= not Din;
   zAdder: IntAdder_8_f400_uid6
      port map ( Cin => cInZ,
                 R => intZout,
                 X => Zin,
                 Y => newAtan2PowStage);

   intDout <= intZout(7);
   Xout <= intXout;
   Yout <= intYout;
   Zout <= intZout;
   Dout <= intDout;
end architecture;

--------------------------------------------------------------------------------
--                           IntAdder_8_f400_uid14
--                      (IntAdderClassical_8_f400_uid16)
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

entity IntAdder_8_f400_uid14 is
   port ( X : in  std_logic_vector(7 downto 0);
          Y : in  std_logic_vector(7 downto 0);
          Cin : in std_logic;
          R : out  std_logic_vector(7 downto 0)   );
end entity;

architecture arch of IntAdder_8_f400_uid14 is
begin
   --Classical
    R <= X + Y + Cin;
end architecture;

--------------------------------------------------------------------------------
--                      FixMicroRotation_8_uid12_stage1
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

entity FixMicroRotation_8_uid12_stage1 is
   port ( Xin : in  std_logic_vector(7 downto 0);
          Yin : in  std_logic_vector(7 downto 0);
          Zin : in  std_logic_vector(7 downto 0);
          Din : in std_logic;
          Xout : out  std_logic_vector(7 downto 0);
          Yout : out  std_logic_vector(7 downto 0);
          Zout : out  std_logic_vector(7 downto 0);
          Dout : out std_logic   );
end entity;

architecture arch of FixMicroRotation_8_uid12_stage1 is
   component IntAdder_8_f400_uid14 is
      port ( X : in  std_logic_vector(7 downto 0);
             Y : in  std_logic_vector(7 downto 0);
             Cin : in std_logic;
             R : out  std_logic_vector(7 downto 0)   );
   end component;

signal Xin_signExtend : std_logic;
signal Xin_shifted :  std_logic_vector(6 downto 0);
signal XinShift :  std_logic_vector(7 downto 0);
signal Yin_signExtend : std_logic;
signal Yin_shifted :  std_logic_vector(6 downto 0);
signal YinShift :  std_logic_vector(7 downto 0);
signal newXinShift :  std_logic_vector(7 downto 0);
signal newYinShift :  std_logic_vector(7 downto 0);
signal cInNewX : std_logic;
signal cInNewY : std_logic;
signal intXout :  std_logic_vector(7 downto 0);
signal intYout :  std_logic_vector(7 downto 0);
signal atan2PowStage :  std_logic_vector(7 downto 0);
signal newAtan2PowStage :  std_logic_vector(7 downto 0);
signal cInZ : std_logic;
signal intZout :  std_logic_vector(7 downto 0);
signal intDout : std_logic;
begin
   Xin_signExtend <= Xin(7);
   Xin_shifted <= Xin(7 downto 1);
   XinShift <= Xin_signExtend & Xin_shifted;
   Yin_signExtend <= Yin(7);
   Yin_shifted <= Yin(7 downto 1);
   YinShift <= Yin_signExtend & Yin_shifted;
   newXinShift <= XinShift xor (7 downto 0 => Din);
   newYinShift <= YinShift xor (7 downto 0 => (not Din));
   cInNewX<= Din;
   cInNewY<= not Din;
   xAdder: IntAdder_8_f400_uid14
      port map ( Cin => cInNewY,
                 R => intXout,
                 X => Xin,
                 Y => newYinShift);

   yAdder: IntAdder_8_f400_uid14
      port map ( Cin => cInNewX,
                 R => intYout,
                 X => Yin,
                 Y => newXinShift);

   atan2PowStage <= '0' & "0000111";
   newAtan2PowStage <= atan2PowStage xor (7 downto 0 => (not Din));
   cInZ<= not Din;
   zAdder: IntAdder_8_f400_uid14
      port map ( Cin => cInZ,
                 R => intZout,
                 X => Zin,
                 Y => newAtan2PowStage);

   intDout <= intZout(7);
   Xout <= intXout;
   Yout <= intYout;
   Zout <= intZout;
   Dout <= intDout;
end architecture;

--------------------------------------------------------------------------------
--                           IntAdder_8_f400_uid22
--                      (IntAdderClassical_8_f400_uid24)
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

entity IntAdder_8_f400_uid22 is
   port ( X : in  std_logic_vector(7 downto 0);
          Y : in  std_logic_vector(7 downto 0);
          Cin : in std_logic;
          R : out  std_logic_vector(7 downto 0)   );
end entity;

architecture arch of IntAdder_8_f400_uid22 is
begin
   --Classical
    R <= X + Y + Cin;
end architecture;

--------------------------------------------------------------------------------
--                      FixMicroRotation_8_uid20_stage2
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

entity FixMicroRotation_8_uid20_stage2 is
   port ( Xin : in  std_logic_vector(7 downto 0);
          Yin : in  std_logic_vector(7 downto 0);
          Zin : in  std_logic_vector(7 downto 0);
          Din : in std_logic;
          Xout : out  std_logic_vector(7 downto 0);
          Yout : out  std_logic_vector(7 downto 0);
          Zout : out  std_logic_vector(7 downto 0);
          Dout : out std_logic   );
end entity;

architecture arch of FixMicroRotation_8_uid20_stage2 is
   component IntAdder_8_f400_uid22 is
      port ( X : in  std_logic_vector(7 downto 0);
             Y : in  std_logic_vector(7 downto 0);
             Cin : in std_logic;
             R : out  std_logic_vector(7 downto 0)   );
   end component;

signal Xin_signExtend :  std_logic_vector(1 downto 0);
signal Xin_shifted :  std_logic_vector(5 downto 0);
signal XinShift :  std_logic_vector(7 downto 0);
signal Yin_signExtend :  std_logic_vector(1 downto 0);
signal Yin_shifted :  std_logic_vector(5 downto 0);
signal YinShift :  std_logic_vector(7 downto 0);
signal newXinShift :  std_logic_vector(7 downto 0);
signal newYinShift :  std_logic_vector(7 downto 0);
signal cInNewX : std_logic;
signal cInNewY : std_logic;
signal intXout :  std_logic_vector(7 downto 0);
signal intYout :  std_logic_vector(7 downto 0);
signal atan2PowStage :  std_logic_vector(7 downto 0);
signal newAtan2PowStage :  std_logic_vector(7 downto 0);
signal cInZ : std_logic;
signal intZout :  std_logic_vector(7 downto 0);
signal intDout : std_logic;
begin
   Xin_signExtend <= (others => Xin(7));
   Xin_shifted <= Xin(7 downto 2);
   XinShift <= Xin_signExtend & Xin_shifted;
   Yin_signExtend <= (others => Yin(7));
   Yin_shifted <= Yin(7 downto 2);
   YinShift <= Yin_signExtend & Yin_shifted;
   newXinShift <= XinShift xor (7 downto 0 => Din);
   newYinShift <= YinShift xor (7 downto 0 => (not Din));
   cInNewX<= Din;
   cInNewY<= not Din;
   xAdder: IntAdder_8_f400_uid22
      port map ( Cin => cInNewY,
                 R => intXout,
                 X => Xin,
                 Y => newYinShift);

   yAdder: IntAdder_8_f400_uid22
      port map ( Cin => cInNewX,
                 R => intYout,
                 X => Yin,
                 Y => newXinShift);

   atan2PowStage <= '0' & "0000100";
   newAtan2PowStage <= atan2PowStage xor (7 downto 0 => (not Din));
   cInZ<= not Din;
   zAdder: IntAdder_8_f400_uid22
      port map ( Cin => cInZ,
                 R => intZout,
                 X => Zin,
                 Y => newAtan2PowStage);

   intDout <= intZout(7);
   Xout <= intXout;
   Yout <= intYout;
   Zout <= intZout;
   Dout <= intDout;
end architecture;

--------------------------------------------------------------------------------
--                           IntAdder_8_f400_uid30
--                      (IntAdderClassical_8_f400_uid32)
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

entity IntAdder_8_f400_uid30 is
   port ( X : in  std_logic_vector(7 downto 0);
          Y : in  std_logic_vector(7 downto 0);
          Cin : in std_logic;
          R : out  std_logic_vector(7 downto 0)   );
end entity;

architecture arch of IntAdder_8_f400_uid30 is
begin
   --Classical
    R <= X + Y + Cin;
end architecture;

--------------------------------------------------------------------------------
--                      FixMicroRotation_8_uid28_stage3
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

entity FixMicroRotation_8_uid28_stage3 is
   port ( Xin : in  std_logic_vector(7 downto 0);
          Yin : in  std_logic_vector(7 downto 0);
          Zin : in  std_logic_vector(7 downto 0);
          Din : in std_logic;
          Xout : out  std_logic_vector(7 downto 0);
          Yout : out  std_logic_vector(7 downto 0);
          Zout : out  std_logic_vector(7 downto 0);
          Dout : out std_logic   );
end entity;

architecture arch of FixMicroRotation_8_uid28_stage3 is
   component IntAdder_8_f400_uid30 is
      port ( X : in  std_logic_vector(7 downto 0);
             Y : in  std_logic_vector(7 downto 0);
             Cin : in std_logic;
             R : out  std_logic_vector(7 downto 0)   );
   end component;

signal Xin_signExtend :  std_logic_vector(2 downto 0);
signal Xin_shifted :  std_logic_vector(4 downto 0);
signal XinShift :  std_logic_vector(7 downto 0);
signal Yin_signExtend :  std_logic_vector(2 downto 0);
signal Yin_shifted :  std_logic_vector(4 downto 0);
signal YinShift :  std_logic_vector(7 downto 0);
signal newXinShift :  std_logic_vector(7 downto 0);
signal newYinShift :  std_logic_vector(7 downto 0);
signal cInNewX : std_logic;
signal cInNewY : std_logic;
signal intXout :  std_logic_vector(7 downto 0);
signal intYout :  std_logic_vector(7 downto 0);
signal atan2PowStage :  std_logic_vector(7 downto 0);
signal newAtan2PowStage :  std_logic_vector(7 downto 0);
signal cInZ : std_logic;
signal intZout :  std_logic_vector(7 downto 0);
signal intDout : std_logic;
begin
   Xin_signExtend <= (others => Xin(7));
   Xin_shifted <= Xin(7 downto 3);
   XinShift <= Xin_signExtend & Xin_shifted;
   Yin_signExtend <= (others => Yin(7));
   Yin_shifted <= Yin(7 downto 3);
   YinShift <= Yin_signExtend & Yin_shifted;
   newXinShift <= XinShift xor (7 downto 0 => Din);
   newYinShift <= YinShift xor (7 downto 0 => (not Din));
   cInNewX<= Din;
   cInNewY<= not Din;
   xAdder: IntAdder_8_f400_uid30
      port map ( Cin => cInNewY,
                 R => intXout,
                 X => Xin,
                 Y => newYinShift);

   yAdder: IntAdder_8_f400_uid30
      port map ( Cin => cInNewX,
                 R => intYout,
                 X => Yin,
                 Y => newXinShift);

   atan2PowStage <= '0' & "0000010";
   newAtan2PowStage <= atan2PowStage xor (7 downto 0 => (not Din));
   cInZ<= not Din;
   zAdder: IntAdder_8_f400_uid30
      port map ( Cin => cInZ,
                 R => intZout,
                 X => Zin,
                 Y => newAtan2PowStage);

   intDout <= intZout(7);
   Xout <= intXout;
   Yout <= intYout;
   Zout <= intZout;
   Dout <= intDout;
end architecture;

--------------------------------------------------------------------------------
--                           IntAdder_8_f400_uid38
--                      (IntAdderClassical_8_f400_uid40)
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

entity IntAdder_8_f400_uid38 is
   port ( X : in  std_logic_vector(7 downto 0);
          Y : in  std_logic_vector(7 downto 0);
          Cin : in std_logic;
          R : out  std_logic_vector(7 downto 0)   );
end entity;

architecture arch of IntAdder_8_f400_uid38 is
begin
   --Classical
    R <= X + Y + Cin;
end architecture;

--------------------------------------------------------------------------------
--                      FixMicroRotation_8_uid36_stage4
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

entity FixMicroRotation_8_uid36_stage4 is
   port ( Xin : in  std_logic_vector(7 downto 0);
          Yin : in  std_logic_vector(7 downto 0);
          Zin : in  std_logic_vector(7 downto 0);
          Din : in std_logic;
          Xout : out  std_logic_vector(7 downto 0);
          Yout : out  std_logic_vector(7 downto 0);
          Zout : out  std_logic_vector(7 downto 0);
          Dout : out std_logic   );
end entity;

architecture arch of FixMicroRotation_8_uid36_stage4 is
   component IntAdder_8_f400_uid38 is
      port ( X : in  std_logic_vector(7 downto 0);
             Y : in  std_logic_vector(7 downto 0);
             Cin : in std_logic;
             R : out  std_logic_vector(7 downto 0)   );
   end component;

signal Xin_signExtend :  std_logic_vector(3 downto 0);
signal Xin_shifted :  std_logic_vector(3 downto 0);
signal XinShift :  std_logic_vector(7 downto 0);
signal Yin_signExtend :  std_logic_vector(3 downto 0);
signal Yin_shifted :  std_logic_vector(3 downto 0);
signal YinShift :  std_logic_vector(7 downto 0);
signal newXinShift :  std_logic_vector(7 downto 0);
signal newYinShift :  std_logic_vector(7 downto 0);
signal cInNewX : std_logic;
signal cInNewY : std_logic;
signal intXout :  std_logic_vector(7 downto 0);
signal intYout :  std_logic_vector(7 downto 0);
signal atan2PowStage :  std_logic_vector(7 downto 0);
signal newAtan2PowStage :  std_logic_vector(7 downto 0);
signal cInZ : std_logic;
signal intZout :  std_logic_vector(7 downto 0);
signal intDout : std_logic;
begin
   Xin_signExtend <= (others => Xin(7));
   Xin_shifted <= Xin(7 downto 4);
   XinShift <= Xin_signExtend & Xin_shifted;
   Yin_signExtend <= (others => Yin(7));
   Yin_shifted <= Yin(7 downto 4);
   YinShift <= Yin_signExtend & Yin_shifted;
   newXinShift <= XinShift xor (7 downto 0 => Din);
   newYinShift <= YinShift xor (7 downto 0 => (not Din));
   cInNewX<= Din;
   cInNewY<= not Din;
   xAdder: IntAdder_8_f400_uid38
      port map ( Cin => cInNewY,
                 R => intXout,
                 X => Xin,
                 Y => newYinShift);

   yAdder: IntAdder_8_f400_uid38
      port map ( Cin => cInNewX,
                 R => intYout,
                 X => Yin,
                 Y => newXinShift);

   atan2PowStage <= '0' & "0000001";
   newAtan2PowStage <= atan2PowStage xor (7 downto 0 => (not Din));
   cInZ<= not Din;
   zAdder: IntAdder_8_f400_uid38
      port map ( Cin => cInZ,
                 R => intZout,
                 X => Zin,
                 Y => newAtan2PowStage);

   intDout <= intZout(7);
   Xout <= intXout;
   Yout <= intYout;
   Zout <= intZout;
   Dout <= intDout;
end architecture;

--------------------------------------------------------------------------------
--                           IntAdder_8_f400_uid46
--                      (IntAdderClassical_8_f400_uid48)
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

entity IntAdder_8_f400_uid46 is
   port ( X : in  std_logic_vector(7 downto 0);
          Y : in  std_logic_vector(7 downto 0);
          Cin : in std_logic;
          R : out  std_logic_vector(7 downto 0)   );
end entity;

architecture arch of IntAdder_8_f400_uid46 is
begin
   --Classical
    R <= X + Y + Cin;
end architecture;

--------------------------------------------------------------------------------
--                      FixMicroRotation_8_uid44_stage5
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

entity FixMicroRotation_8_uid44_stage5 is
   port ( Xin : in  std_logic_vector(7 downto 0);
          Yin : in  std_logic_vector(7 downto 0);
          Zin : in  std_logic_vector(7 downto 0);
          Din : in std_logic;
          Xout : out  std_logic_vector(7 downto 0);
          Yout : out  std_logic_vector(7 downto 0);
          Zout : out  std_logic_vector(7 downto 0);
          Dout : out std_logic   );
end entity;

architecture arch of FixMicroRotation_8_uid44_stage5 is
   component IntAdder_8_f400_uid46 is
      port ( X : in  std_logic_vector(7 downto 0);
             Y : in  std_logic_vector(7 downto 0);
             Cin : in std_logic;
             R : out  std_logic_vector(7 downto 0)   );
   end component;

signal Xin_signExtend :  std_logic_vector(4 downto 0);
signal Xin_shifted :  std_logic_vector(2 downto 0);
signal XinShift :  std_logic_vector(7 downto 0);
signal Yin_signExtend :  std_logic_vector(4 downto 0);
signal Yin_shifted :  std_logic_vector(2 downto 0);
signal YinShift :  std_logic_vector(7 downto 0);
signal newXinShift :  std_logic_vector(7 downto 0);
signal newYinShift :  std_logic_vector(7 downto 0);
signal cInNewX : std_logic;
signal cInNewY : std_logic;
signal intXout :  std_logic_vector(7 downto 0);
signal intYout :  std_logic_vector(7 downto 0);
signal atan2PowStage :  std_logic_vector(7 downto 0);
signal newAtan2PowStage :  std_logic_vector(7 downto 0);
signal cInZ : std_logic;
signal intZout :  std_logic_vector(7 downto 0);
signal intDout : std_logic;
begin
   Xin_signExtend <= (others => Xin(7));
   Xin_shifted <= Xin(7 downto 5);
   XinShift <= Xin_signExtend & Xin_shifted;
   Yin_signExtend <= (others => Yin(7));
   Yin_shifted <= Yin(7 downto 5);
   YinShift <= Yin_signExtend & Yin_shifted;
   newXinShift <= XinShift xor (7 downto 0 => Din);
   newYinShift <= YinShift xor (7 downto 0 => (not Din));
   cInNewX<= Din;
   cInNewY<= not Din;
   xAdder: IntAdder_8_f400_uid46
      port map ( Cin => cInNewY,
                 R => intXout,
                 X => Xin,
                 Y => newYinShift);

   yAdder: IntAdder_8_f400_uid46
      port map ( Cin => cInNewX,
                 R => intYout,
                 X => Yin,
                 Y => newXinShift);

   atan2PowStage <= '0' & "0000000";
   newAtan2PowStage <= atan2PowStage xor (7 downto 0 => (not Din));
   cInZ<= not Din;
   zAdder: IntAdder_8_f400_uid46
      port map ( Cin => cInZ,
                 R => intZout,
                 X => Zin,
                 Y => newAtan2PowStage);

   intDout <= intZout(7);
   Xout <= intXout;
   Yout <= intYout;
   Zout <= intZout;
   Dout <= intDout;
end architecture;

--------------------------------------------------------------------------------
--     FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_0
-- This operator is part of the Infinite Virtual Library FloPoCoLib
-- All rights reserved 
-- Authors: Florent de Dinechin (2007)
--------------------------------------------------------------------------------
library ieee; 
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
library work;
entity FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_0 is
   port ( X : in  std_logic_vector(3 downto 0);
          Y : out  std_logic_vector(3 downto 0)   );
end entity;

architecture arch of FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_0 is
begin
  with X select  Y <= 
   "0000" when "0000",
   "0001" when "0001",
   "0001" when "0010",
   "0010" when "0011",
   "0010" when "0100",
   "0011" when "0101",
   "0100" when "0110",
   "0100" when "0111",
   "0101" when "1000",
   "0101" when "1001",
   "0110" when "1010",
   "0111" when "1011",
   "0111" when "1100",
   "1000" when "1101",
   "1001" when "1110",
   "1001" when "1111",
   "----" when others;
end architecture;

--------------------------------------------------------------------------------
--     FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_1
-- This operator is part of the Infinite Virtual Library FloPoCoLib
-- All rights reserved 
-- Authors: Florent de Dinechin (2007)
--------------------------------------------------------------------------------
library ieee; 
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
library work;
entity FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_1 is
   port ( X : in  std_logic_vector(4 downto 0);
          Y : out  std_logic_vector(8 downto 0)   );
end entity;

architecture arch of FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_1 is
begin
  with X select  Y <= 
   "000000000" when "00000",
   "000001010" when "00001",
   "000010011" when "00010",
   "000011101" when "00011",
   "000100111" when "00100",
   "000110001" when "00101",
   "000111010" when "00110",
   "001000100" when "00111",
   "001001110" when "01000",
   "001010111" when "01001",
   "001100001" when "01010",
   "001101011" when "01011",
   "001110101" when "01100",
   "001111110" when "01101",
   "010001000" when "01110",
   "010010010" when "01111",
   "101100101" when "10000",
   "101101110" when "10001",
   "101111000" when "10010",
   "110000010" when "10011",
   "110001011" when "10100",
   "110010101" when "10101",
   "110011111" when "10110",
   "110101001" when "10111",
   "110110010" when "11000",
   "110111100" when "11001",
   "111000110" when "11010",
   "111001111" when "11011",
   "111011001" when "11100",
   "111100011" when "11101",
   "111101101" when "11110",
   "111110110" when "11111",
   "---------" when others;
end architecture;

--------------------------------------------------------------------------------
--                           IntAdder_9_f400_uid55
--                      (IntAdderClassical_9_f400_uid57)
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

entity IntAdder_9_f400_uid55 is
   port ( X : in  std_logic_vector(8 downto 0);
          Y : in  std_logic_vector(8 downto 0);
          Cin : in std_logic;
          R : out  std_logic_vector(8 downto 0)   );
end entity;

architecture arch of IntAdder_9_f400_uid55 is
begin
   --Classical
    R <= X + Y + Cin;
end architecture;

--------------------------------------------------------------------------------
--         FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed
-- This operator is part of the Infinite Virtual Library FloPoCoLib
-- All rights reserved 
-- Authors: 
--------------------------------------------------------------------------------
-- combinatorial

library ieee;
use ieee.std_logic_1164.all;
use ieee.std_logic_arith.all;
use ieee.std_logic_unsigned.all;
library std;
use std.textio.all;
library work;

entity FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed is
   port ( X : in  std_logic_vector(8 downto 0);
          R : out  std_logic_vector(8 downto 0)   );
end entity;

architecture arch of FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed is
   component FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_0 is
      port ( X : in  std_logic_vector(3 downto 0);
             Y : out  std_logic_vector(3 downto 0)   );
   end component;

   component FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_1 is
      port ( X : in  std_logic_vector(4 downto 0);
             Y : out  std_logic_vector(8 downto 0)   );
   end component;

   component IntAdder_9_f400_uid55 is
      port ( X : in  std_logic_vector(8 downto 0);
             Y : in  std_logic_vector(8 downto 0);
             Cin : in std_logic;
             R : out  std_logic_vector(8 downto 0)   );
   end component;

signal d0 :  std_logic_vector(3 downto 0);
signal pp0 :  std_logic_vector(3 downto 0);
signal addOp0 :  std_logic_vector(8 downto 0);
signal d1 :  std_logic_vector(4 downto 0);
signal pp1 :  std_logic_vector(8 downto 0);
signal addOp1 :  std_logic_vector(8 downto 0);
signal OutRes :  std_logic_vector(8 downto 0);
attribute rom_extract: string;
attribute rom_style: string;
attribute rom_extract of FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_0: component is "yes";
attribute rom_extract of FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_1: component is "yes";
attribute rom_style of FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_0: component is "distributed";
attribute rom_style of FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_1: component is "distributed";
begin
   d0 <= X(3 downto 0);
   KCMTable_0: FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_0
      port map ( X => d0,
                 Y => pp0);
   addOp0 <= (8 downto 4 => '0') & pp0;
   d1 <= X(8 downto 4);
   KCMTable_1: FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed_Table_1
      port map ( X => d1,
                 Y => pp1);
   addOp1 <= pp1;
   Result_Adder: IntAdder_9_f400_uid55
      port map ( Cin => '0',
                 R => OutRes,
                 X => addOp0,
                 Y => addOp1);
   R <= OutRes(8 downto 0);
end architecture;

--------------------------------------------------------------------------------
--                            CordicSinCos_8_uid2
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

entity CordicSinCos_8_uid2 is
   port ( Z : in  std_logic_vector(5 downto 0);
          Xout : out  std_logic_vector(5 downto 0);
          Yout : out  std_logic_vector(5 downto 0)   );
end entity;

architecture arch of CordicSinCos_8_uid2 is
   component FixMicroRotation_8_uid12_stage1 is
      port ( Xin : in  std_logic_vector(7 downto 0);
             Yin : in  std_logic_vector(7 downto 0);
             Zin : in  std_logic_vector(7 downto 0);
             Din : in std_logic;
             Xout : out  std_logic_vector(7 downto 0);
             Yout : out  std_logic_vector(7 downto 0);
             Zout : out  std_logic_vector(7 downto 0);
             Dout : out std_logic   );
   end component;

   component FixMicroRotation_8_uid20_stage2 is
      port ( Xin : in  std_logic_vector(7 downto 0);
             Yin : in  std_logic_vector(7 downto 0);
             Zin : in  std_logic_vector(7 downto 0);
             Din : in std_logic;
             Xout : out  std_logic_vector(7 downto 0);
             Yout : out  std_logic_vector(7 downto 0);
             Zout : out  std_logic_vector(7 downto 0);
             Dout : out std_logic   );
   end component;

   component FixMicroRotation_8_uid28_stage3 is
      port ( Xin : in  std_logic_vector(7 downto 0);
             Yin : in  std_logic_vector(7 downto 0);
             Zin : in  std_logic_vector(7 downto 0);
             Din : in std_logic;
             Xout : out  std_logic_vector(7 downto 0);
             Yout : out  std_logic_vector(7 downto 0);
             Zout : out  std_logic_vector(7 downto 0);
             Dout : out std_logic   );
   end component;

   component FixMicroRotation_8_uid36_stage4 is
      port ( Xin : in  std_logic_vector(7 downto 0);
             Yin : in  std_logic_vector(7 downto 0);
             Zin : in  std_logic_vector(7 downto 0);
             Din : in std_logic;
             Xout : out  std_logic_vector(7 downto 0);
             Yout : out  std_logic_vector(7 downto 0);
             Zout : out  std_logic_vector(7 downto 0);
             Dout : out std_logic   );
   end component;

   component FixMicroRotation_8_uid44_stage5 is
      port ( Xin : in  std_logic_vector(7 downto 0);
             Yin : in  std_logic_vector(7 downto 0);
             Zin : in  std_logic_vector(7 downto 0);
             Din : in std_logic;
             Xout : out  std_logic_vector(7 downto 0);
             Yout : out  std_logic_vector(7 downto 0);
             Zout : out  std_logic_vector(7 downto 0);
             Dout : out std_logic   );
   end component;

   component FixMicroRotation_8_uid4_stage0 is
      port ( Xin : in  std_logic_vector(7 downto 0);
             Yin : in  std_logic_vector(7 downto 0);
             Zin : in  std_logic_vector(7 downto 0);
             Din : in std_logic;
             Xout : out  std_logic_vector(7 downto 0);
             Yout : out  std_logic_vector(7 downto 0);
             Zout : out  std_logic_vector(7 downto 0);
             Dout : out std_logic   );
   end component;

   component FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed is
      port ( X : in  std_logic_vector(8 downto 0);
             R : out  std_logic_vector(8 downto 0)   );
   end component;

signal X0 :  std_logic_vector(7 downto 0);
signal Y0 :  std_logic_vector(7 downto 0);
signal Z0 :  std_logic_vector(7 downto 0);
signal D0 : std_logic;
signal X1 :  std_logic_vector(7 downto 0);
signal Y1 :  std_logic_vector(7 downto 0);
signal Z1 :  std_logic_vector(7 downto 0);
signal D1 : std_logic;
signal X2 :  std_logic_vector(7 downto 0);
signal Y2 :  std_logic_vector(7 downto 0);
signal Z2 :  std_logic_vector(7 downto 0);
signal D2 : std_logic;
signal X3 :  std_logic_vector(7 downto 0);
signal Y3 :  std_logic_vector(7 downto 0);
signal Z3 :  std_logic_vector(7 downto 0);
signal D3 : std_logic;
signal X4 :  std_logic_vector(7 downto 0);
signal Y4 :  std_logic_vector(7 downto 0);
signal Z4 :  std_logic_vector(7 downto 0);
signal D4 : std_logic;
signal X5 :  std_logic_vector(7 downto 0);
signal Y5 :  std_logic_vector(7 downto 0);
signal Z5 :  std_logic_vector(7 downto 0);
signal D5 : std_logic;
signal X6 :  std_logic_vector(7 downto 0);
signal Y6 :  std_logic_vector(7 downto 0);
signal Z6 :  std_logic_vector(7 downto 0);
signal D6 : std_logic;
signal preXout :  std_logic_vector(8 downto 0);
signal intXout :  std_logic_vector(8 downto 0);
signal preYout :  std_logic_vector(8 downto 0);
signal intYout :  std_logic_vector(8 downto 0);
signal roundedIntXout :  std_logic_vector(6 downto 0);
signal roundedIntYout :  std_logic_vector(6 downto 0);
begin
   X0<= '0' & "0010000";
   Y0<= '0' & "0000000";
   Z0<= Z & "00";
   D0<= Z(5);
   microRotation0: FixMicroRotation_8_uid4_stage0
      port map ( Din => D0,
                 Dout => D1,
                 Xin => X0,
                 Xout => X1,
                 Yin => Y0,
                 Yout => Y1,
                 Zin => Z0,
                 Zout => Z1);

   microRotation1: FixMicroRotation_8_uid12_stage1
      port map ( Din => D1,
                 Dout => D2,
                 Xin => X1,
                 Xout => X2,
                 Yin => Y1,
                 Yout => Y2,
                 Zin => Z1,
                 Zout => Z2);

   microRotation2: FixMicroRotation_8_uid20_stage2
      port map ( Din => D2,
                 Dout => D3,
                 Xin => X2,
                 Xout => X3,
                 Yin => Y2,
                 Yout => Y3,
                 Zin => Z2,
                 Zout => Z3);

   microRotation3: FixMicroRotation_8_uid28_stage3
      port map ( Din => D3,
                 Dout => D4,
                 Xin => X3,
                 Xout => X4,
                 Yin => Y3,
                 Yout => Y4,
                 Zin => Z3,
                 Zout => Z4);

   microRotation4: FixMicroRotation_8_uid36_stage4
      port map ( Din => D4,
                 Dout => D5,
                 Xin => X4,
                 Xout => X5,
                 Yin => Y4,
                 Yout => Y5,
                 Zin => Z4,
                 Zout => Z5);

   microRotation5: FixMicroRotation_8_uid44_stage5
      port map ( Din => D5,
                 Dout => D6,
                 Xin => X5,
                 Xout => X6,
                 Yin => Y5,
                 Yout => Y6,
                 Zin => Z5,
                 Zout => Z6);

   preXout<= X6 & '0';
   constMultiplierX: FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed
      port map ( R => intXout,
                 X => preXout);

   preYout<= Y6 & '0';
   constMultiplierY: FixRealKCM_M4_3_M4_0_607252935008881256169446834473_signed
      port map ( R => intYout,
                 X => preYout);

   roundedIntXout<= intXout(8 downto 2) + ("000000" & '1');
   roundedIntYout<= intYout(8 downto 2) + ("000000" & '1');
   Xout<= roundedIntXout(6 downto 1);
   Yout<= roundedIntYout(6 downto 1);
end architecture;

--------------------------------------------------------------------------------
--                       TestBench_CordicSinCos_8_uid2
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

entity TestBench_CordicSinCos_8_uid2 is
end entity;

architecture behavorial of TestBench_CordicSinCos_8_uid2 is
   component CordicSinCos_8_uid2 is
      port ( Z : in  std_logic_vector(5 downto 0);
             Xout : out  std_logic_vector(5 downto 0);
             Yout : out  std_logic_vector(5 downto 0)   );
   end component;
   signal Z :  std_logic_vector(5 downto 0);
   signal Xout :  std_logic_vector(5 downto 0);
   signal Yout :  std_logic_vector(5 downto 0);
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
   test: CordicSinCos_8_uid2
      port map ( Xout => Xout,
                 Yout => Yout,
                 Z => Z);
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
      variable V_Z : bit_vector(5 downto 0);
      variable V_Xout : bit_vector(5 downto 0);
      variable V_Yout : bit_vector(5 downto 0);
   begin
      -- Send reset
      rst <= '1';
      wait for 10 ns;
      rst <= '0';
      while not endfile(inputsFile) loop
          -- positionning inputs
         readline(inputsFile,inline);
         read(inline ,V_Z);
         read(inline,tmpChar);
         Z <= to_stdlogicvector(V_Z);
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
      variable V_Z : bit_vector(5 downto 0);
      variable V_Xout : bit_vector(5 downto 0);
      variable V_Yout : bit_vector(5 downto 0);
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
          wait for 10 ns; -- wait for pipeline to flush
         counter := counter + 2;
      end loop;
      report (integer'image(errorCounter) & " error(s) encoutered.");
      assert false report "End of simulation" severity failure;
   end process;

end architecture;

