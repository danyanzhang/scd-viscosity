function lines = drawselection()
% A figure must be open for this to work

i = 1; % counter
collectInput = true;

uiwait(msgbox('Draw a line to select an individual oxygen cycle'))
while collectInput == true;
    lines{i} = drawline;
    userInput = questdlg('Continue data selection?', 'Attention', 'New Line', 'Finish and Save', 'New Line');
    switch userInput
    case 'New Line'
        i = i+1;
    case 'Finish and Save'
        collectInput = false;
    end
end

end