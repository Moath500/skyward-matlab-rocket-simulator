function data = RocketGeometryGUI
% Rocket and Fins geometry GUI to develop DATCOM file
%
% How to use: 
% Run the script and simply modify the default values.
% The plots are updated automatically once the numerical values are changed
% Press the button to printf the DATCOM geomtry values in console
%
% Author: Luca Facchini
% Skyward Experimental Rocketry | MSA Dept
% email: luca.facchini@skywarder.eu
% Website: http://www.skywarder.eu
% Release date: 11 October 2019 | First Version
% License:  2-clause BSD

close all
clear all
clc
%% Set default values - Dimensions in cm
data = struct();
data.D = 9;
data.RocketLength = 190;
data.NoseLength = 30;
data.BodyLength = data.RocketLength - data.NoseLength;
data.FinMaxChord = 10;
data.FinMinChord = 5;
data.FinHeight = 5;
data.BottomDist = 0.0;
data.XLe = data.RocketLength - data.BottomDist - [data.FinMaxChord ; ...
    data.FinMaxChord-(data.FinMaxChord-data.FinMinChord)/2];
data.XCG = 126;

% Fin cross section
data.FinT = 0.4;
data.LmaxMaxChord = 0.6;
data.LmaxMinChord = 0.6;
data.LflatMaxChord = 8;
data.LflatMinChord = 4;
data.LmaxMaxChord = 0.5*(data.FinMaxChord - data.LflatMaxChord);
data.LmaxMinChord = 0.5*(data.FinMinChord - data.LflatMinChord);

%% Create layout
f = figure('Visible','off','Position',[100 100 763 638]);
f.Name = 'Rocket Geometry';

RocketUIAxes = axes('Units','Pixels','Position',[56 355 503 265],'Title',...
    'Rocket');
RocketUIAxes.XGrid = 'on';
RocketUIAxes.YGrid = 'on';
xlabel(RocketUIAxes, 'X [cm]')
ylabel(RocketUIAxes, 'Y [cm]')
title(RocketUIAxes,'Rocket');

FinsUIAxes = axes('Units','Pixels','Position',[56 50 503 255],'Title',...
    'Fins');
FinsUIAxes.XGrid = 'on';
FinsUIAxes.YGrid = 'on';
xlabel(FinsUIAxes, 'X [cm]')
ylabel(FinsUIAxes, 'Y [cm]')
title(FinsUIAxes,'Fins');


Label = uicontrol('Style','text','String','Measures in [cm]',...
    'Position',[636 595 100 22]);

RocketLengthLabel = uicontrol('Style','text','String','Rocket Length',...
    'Position',[575 566 83 22],'Parent',f);
RocketLengthEditField = uicontrol('Style','edit',...
    'Position',[660 570 83 22 ]);

NoseLengthLabel = uicontrol('Style','text','String','Nose Length',...
    'Position',[575 536 83 22]);
NoseLengthEditField = uicontrol('Style','edit',...
    'Position',[660 540 83 22 ]);

MaxChordLabel = uicontrol('Style','text','String','Max Chord',...
    'Position',[575 506 83 22]);
MaxChordEditField = uicontrol('Style','edit',...
    'Position',[660 510 83 22 ]);

MinChordLabel = uicontrol('Style','text','String','Max Chord',...
    'Position',[575 476 83 22]);
MinChordEditField = uicontrol('Style','edit',...
    'Position',[660 480 83 22 ]);

XLe1Label = uicontrol('Style','text','String','XLe(1)',...
    'Position',[575 446 83 22]);
XLe1EditField = uicontrol('Style','edit',...
    'Position',[660 450 83 22 ]);

XLe2Label = uicontrol('Style','text','String','XLe(2)',...
    'Position',[575 416 83 22]);
XLe2EditField = uicontrol('Style','edit',...
    'Position',[660 420 83 22 ]);

FinHLabel = uicontrol('Style','text','String','Fin Height',...
    'Position',[575 386 83 22]);
FinHEditField = uicontrol('Style','edit',...
    'Position',[660 390 83 22 ]);

XCGLabel = uicontrol('Style','text','String','XCG',...
    'Position',[575 356 83 22]);
XCGEditField = uicontrol('Style','edit',...
    'Position',[660 360 83 22 ]);


ThicknessLabel = uicontrol('Style','text','String','Thickness',...
    'Position',[562 246 83 22]);
ThicknessEditField = uicontrol('Style','edit',...
    'Position',[660 250 83 22 ]);

LflatMaxCLabel = uicontrol('Style','text','String','Lflat (Max chord)',...
    'Position',[562 216 103 22]);
LflatMaxCEditField = uicontrol('Style','edit',...
    'Position',[660 220 83 22 ]);
LflatMaxCEditField.Callback = @LflatMaxCEditFieldValueChanged;

LflatMinCLabel = uicontrol('Style','text','String','Lflat (Min chord)',...
    'Position',[562 186 103 22]);
LflatMinCEditField = uicontrol('Style','edit',...
    'Position',[660 190 83 22 ]);

Button = uicontrol('Style','pushbutton','String','Print DATCOM text',...
    'Position',[630 145 100 22 ],'Callback',@PrintDATCOM);


%% Set default values
RocketLengthEditField.String = num2str(data.RocketLength);
NoseLengthEditField.String = num2str(data.NoseLength);
MaxChordEditField.String = num2str(data.FinMaxChord);
MinChordEditField.String = num2str(data.FinMinChord);
XLe1EditField.String = num2str(data.XLe(1));
XLe2EditField.String = num2str(data.XLe(2));
FinHEditField.String = num2str(data.FinHeight);
XCGEditField.String = num2str(data.XCG);

ThicknessEditField.String = num2str(data.FinT);
LflatMaxCEditField.String = num2str(data.LflatMaxChord);
LflatMinCEditField.String = num2str(data.LflatMinChord);

updateplot;

%% Set Callbacks

NoseLengthEditField.Callback = @NoseLengthEditFieldValueChanged;
RocketLengthEditField.Callback = @RocketLengthEditFieldValueChanged;
XCGEditField.Callback = @XCGEditFieldValueChanged;
LflatMaxCEditField.Callback = @LflatMaxCEditFieldValueChanged;
LflatMinCEditField.Callback = @LflatMinCEditFieldValueChanged;
XLe1EditField.Callback = @XLe1EditFieldValueChanged;
XLe2EditField.Callback = @XLe2EditFieldValueChanged;
MaxChordEditField.Callback = @MaxChordEditFieldValueChanged;
MinChordEditField.Callback = @MinChordEditFieldValueChanged;
FinHEditField.Callback = @FinHeightEditFieldValueChanged;
ThicknessEditField.Callback = @FinThicknessEditFieldValueChanged;

% Make the UI visible.
f.Visible = 'on';

%% FUNCTIONS EVENETS
    function PrintDATCOM(src,event)
        % Print DATCOM output
        str = sprintf("$REFQ\n");
        str = sprintf('%s XCG=%4.3f,\n',str,data.XCG/100);
        str = sprintf('%s SREF=%6.4f,\n',str,pi*(data.D/100)^2/4);
        str = sprintf('%s LREF=%4.3f,\n',str,data.D/100);
        str = sprintf('%s LATREF=%4.2f,$\n',str,data.D/100);
        str = sprintf('%s$AXIBOD\n',str);
        str = sprintf('%s TNOSE=KARMAN,\n',str);
        str = sprintf('%s LNOSE=%4.2f,\n',str,data.NoseLength/100);
        str = sprintf('%s DNOSE=%4.2f,\n',str,data.D/100);
        str = sprintf('%s LCENTR=%4.2f,\n',str,(data.RocketLength-data.NoseLength)/100);
        str = sprintf('%s DCENTR=%4.2f,\n',str,data.D/100);
        str = sprintf('%s DEXIT=%4.2f,\n',str,data.D/100);
        str = sprintf('%s BASE=.FALSE.,$\n',str);
        str = sprintf('%s$FINSET1\n',str);
        str = sprintf('%s NPANEL=4.0,\n',str);
        str = sprintf('%s XLE=%4.3f,%4.3f,\n',str,data.XLe(1)/100,data.XLe(2)/100);
        str = sprintf('%s PHIF=0.0,90.0,180.0,270.0,\n',str);
        str = sprintf('%s LER=2*0.0025,\n',str);
        str = sprintf('%s STA=0.0,\n',str);
        str = sprintf('%s SSPAN=%4.3f,%4.3f,\n',str,data.D/2/100,(data.D/2+data.FinHeight)/100);
        str = sprintf('%s CHORD=%4.3f,%4.3f,\n',str,data.FinMaxChord/2/100,data.FinMinChord/100);
        str = sprintf('%s SECTYP=HEX,\n',str);
        str = sprintf('%s ZUPPER=%4.3f,%4.3f,\n',str,data.FinT/2/data.FinMaxChord,data.FinT/2/data.FinMinChord);
        str = sprintf('%s LMAXU=%4.3f,%4.3f,\n',str,data.LmaxMaxChord/data.FinMaxChord,...
            data.LmaxMinChord/data.FinMinChord);
        str = sprintf('%s LFLATU=%3.2f,%3.2f,$\n',str,data.LflatMaxChord/data.FinMaxChord,...
            data.LflatMinChord/data.FinMinChord);
             
        fprintf(str)
    end

     % Value changed function: RocketLengthEditField
    function RocketLengthEditFieldValueChanged(src,event)
        value = str2double(src.String);
        data.RocketLength = value;
        updateplot;
    end

    % Value changed function: NoseLengthEditField
    function NoseLengthEditFieldValueChanged(src,event)
        value = str2double(src.String);
        if value > data.RocketLength
            msgbox('Invalid value of nose length. Cannot be > than rocket length','Warning');
            % Reassign previous value
            src.String = num2str(data.NoseLength);
        else % update value
            data.NoseLength = value;
            updateplot;
        end
     end
    
    % Value changed function: XLe1EditField
    function XLe1EditFieldValueChanged(src, event)
        value = str2double(src.String);
        if value > data.RocketLength(1)
            msgbox('Invalid value of XLE(1). Cannot be bigger than rocket length','Warning');
            src.String = num2str(data.XLe(1));
        elseif value > data.XLe(2)
            msgbox('Invalid value of XLE. Cannot be bigger than rocket length','Warning');
            % Reassign previous value
            src.String = num2str(data.XLe(1));
        else
            data.XLe(1) = value;
            updateplot;
        end
    end

    % Value changed function: XLe2EditField
    function XLe2EditFieldValueChanged(src,event)
        value = str2double(src.String);
        % check if value is ok
        if value < data.XLe(1)
            msgbox('Invalid value of XLe(2). Cannot be smaller than XLe(1)','Warning');
            % Reassign previous value
            src.String = num2str(data.XLe(2));
        else 
            data.XLe(2) = value;
            updateplot;
        end
    end

    % Value changed function: MaxChordEditField
    function MaxChordEditFieldValueChanged(src,event)
        value = str2double(src.String);
        if value > data.RocketLength
            msgbox('Invalid value of chord. Cannot be bigger than rocket length','Warning');
            % Reassign previous value
            src.String = num2str(data.FinMaxChord);
        elseif value <= data.FinMinChord
            msgbox('Invalid value of chord. Cannot be smaller than the minor chord','Warning');
            src.String = num2str(data.FinMaxChord);
        else
            data.FinMaxChord = value;
            data.LmaxMaxChord = (data.FinMaxChord - data.LflatMaxChord)/2;
            updateplot;
        end
    end

    % Value changed function: MinChordEditField
    function MinChordEditFieldValueChanged(src,event)
        value = str2double(src.String);
        
        if value > data.RocketLength
            msgbox('Invalid value of chord. Cannot be > than rocket length','Warning');
            % Reassign previous value
            src.String = num2str(data.FinMinChord); 
        elseif value >= data.FinMaxChord
            msgbox('Invalid value of chord. Cannot be > than the minor chord','Warning');
            % Reassign previous value
            src.String = num2str(data.FinMinChord);
        else
            data.FinMinChord = value;
            data.LmaxMinChord = (data.FinMinChord - data.LflatMinChord)/2;
            updateplot;
        end
    end

    % Value changed function: FinHeightEditField
    function FinHeightEditFieldValueChanged(src,event)
        value = str2double(get(src,'String'));
        data.FinHeight = value;
        updateplot;
    end

    % Value changed function: FinThicknessEditField
    function FinThicknessEditFieldValueChanged(src,event)
        value = str2double(get(src,'String'));
        data.FinT = value;
        updateplot;
    end

    % Value changed function: LflatMaxChordEditField
    function LflatMaxCEditFieldValueChanged(src,event)
        value = str2double(get(src,'String'));
        
        if value > data.FinMaxChord
            msgbox('Invalid values of Lmax. Cannot be > than chord','Warning');
            src.String = num2str(data.LflatMaxChord);
        else
            data.LflatMaxChord = value;
            data.LmaxMaxChord = (data.FinMaxChord-data.LflatMaxChord)/2;
            updateplot;
        end
    end

    % Value changed function: LflatMinChordEditField
    function LflatMinCEditFieldValueChanged(src,event)
        value = str2double(get(src,'String'));
        if value > data.FinMinChord
            msgbox('Invalid values of Lmax. Cannot be bigger than chord','Warning');
            src.String = num2str(data.LflatMinChord);
        else
            data.LflatMinChord = value;
            data.LmaxMinChord = (data.FinMinChord-data.LflatMinChord)/2;
            updateplot;
        end
    end

    % Value changed function: XCGEditField
    function XCGEditFieldValueChanged(src,event)
        value = str2double(get(src,'String'));
        data.XCG = value;
        updateplot;
    end

    function updateplot
        % Nose
        xn = linspace(0,data.NoseLength);
        theta = acos(1-2.*xn./data.NoseLength);
        yn = data.D/2./sqrt(pi).*sqrt(theta-sin(2*theta)/2);
        % Body coordinates
        xb = [xn data.RocketLength data.RocketLength];
        yb = [yn data.D/2 0];
        % Fins coordinates
        xf = [data.XLe(1)  data.XLe(1)+data.FinMaxChord data.XLe(2)+data.FinMinChord ...
            data.XLe(2) data.XLe(1)];
        yf = [data.D/2 data.D/2 data.D/2+data.FinHeight data.D/2+data.FinHeight data.D/2];
        h1 = plot(RocketUIAxes,xb,yb,'-b');
        hold(RocketUIAxes,'on');
        plot(RocketUIAxes,xb,-yb,'-b');
        plot(RocketUIAxes,xf,yf,'-b');
        plot(RocketUIAxes,xf,-yf,'-b');
        h2 = plot(RocketUIAxes,data.XCG,0,'or');
        plot(RocketUIAxes,data.XCG,0,'+r');
        RocketUIAxes.YLim = [-data.D*2 data.D*2];
        RocketUIAxes.XLim = [-15 (15+data.RocketLength)];
        xlabel(RocketUIAxes, 'X [cm]')
        ylabel(RocketUIAxes, 'Y [cm]')
        legend([h1,h2],'Rocket','CG');
        hold(RocketUIAxes,'off');
        
        % Draw fins
        xf = [-data.FinMaxChord/2 , -data.FinMaxChord/2+data.LmaxMaxChord , ...
            data.FinMaxChord/2-data.LmaxMaxChord , data.FinMaxChord/2];
        yf = [0 data.FinT/2 data.FinT/2 0];
        h1 = plot(FinsUIAxes,xf,yf,'-b');
        hold(FinsUIAxes,'on');
        plot(FinsUIAxes,xf,-yf,'b');
        
        xf = [-data.FinMinChord/2 , -data.FinMinChord/2+data.LmaxMinChord , ...
            data.FinMinChord/2-data.LmaxMinChord , data.FinMinChord/2];
        yf = [0 data.FinT/2 data.FinT/2 0];
        h2 = plot(FinsUIAxes,xf,yf,'-r');
        plot(FinsUIAxes,xf,-yf,'r');
        FinsUIAxes.XLim = [-data.FinMaxChord/2-0.1 data.FinMaxChord/2+0.1];
        FinsUIAxes.YLim = [-data.FinT data.FinT];
        xlabel(FinsUIAxes, 'X [cm]')
        ylabel(FinsUIAxes, 'Y [cm]')
         legend([h1,h2],'Cross section at fin root','Cross section at fin tip');
        hold(FinsUIAxes,'off');
    end
end
 
