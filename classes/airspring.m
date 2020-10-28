classdef airspring
   properties
       load
       naturalfrequency
       volume
       maxdiameter
       mass
   end
   properties %(Access=private)
      springid=[0:5];
      stylenumber={'29','153-2','119','233-2','28','22'};
      designheight=[9.5,6,5,11.25,9.5,9.5];
      load_40psig=[9780,4564,5490,3413,4590,2060];
      load_60psig=[14860,7048,8450,5631,7010,3170];
      load_80psig=[20060,9682,11450,7691,9590,4280];
      natural_freq=[92,121,138,89,101,106];
      isolation_400cpm=[94.4,89.9,86.5,95,93.2,92.4];
      isolation_800cpm=[98.7,97.9,96.9,98.8,98.4,98.2];
      isolation_1500cpm=[99.6,99.3,99.1,99.7,99.5,99.5];
      volume_100psig=[2934,981,787,1418,1596,782];
      max_height=[10,7,6,13,10,10];
      min_height=[4,2.5,2,3,3.4,3.2];
      max_diameter=[22.7,18.1,17.4,15.5,17.4,12.9];
      brand={'Firestone','Firestone','Firestone','Firestone','Firestone','Firestone'};
     database
   end
   methods
       function database=get.database(obj)
            database=table(obj.springid',obj.stylenumber',obj.designheight',obj.load_40psig',obj.load_60psig',obj.load_80psig',obj.natural_freq',obj.isolation_400cpm',...
       obj.isolation_800cpm',obj.isolation_1500cpm',obj.volume_100psig',obj.max_height',obj.min_height',obj.max_diameter',obj.brand'); 
       end
   end
end