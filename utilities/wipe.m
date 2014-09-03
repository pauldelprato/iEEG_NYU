function wipe

%take user input to make sure you want to (y/n)...
result = input('You sure you want to wipe? (y/n)','s');

if strcmpi(result, 'y')
    fclose all; close all; clear all;
else
    display('Not wiping.')
end