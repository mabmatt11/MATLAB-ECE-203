function display(f)
%DISPLAY Displays the PULSE objects data.
%   DISPLAY(f) displays the objects data.  It is automatically called from
%   when the object is entered in the command window without terminating semicolon.
%
%   See also PULSE

% Jordan Rosenthal, 12/16/97

Param_Names = char( {'Amplitude: ', 'Length: ', 'Delay: '} );
Param_Vals = [f.Amplitude, f.Length, f.Delay];
Divider(1:length(f.Name)) = '-';
disp(' ');
disp(f.Name);
disp(Divider);
for i = 1:size(Param_Names,1)
   disp([Param_Names(i,:) num2str(Param_Vals(i))]);
end
disp(' ');
