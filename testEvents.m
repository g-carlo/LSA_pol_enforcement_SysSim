clear all
close all
clear classes;

asd = eventSLS_Cls('ciao',2,'hello');


field1 = 'type';  
value1 = cell(1,2);
field2 = 'time';  
value2 = cell(1,2);


for n = 1:3
   eventQueue(n) =  eventSLS_Cls('ciao',2,'hello');
end

