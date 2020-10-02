function varargout = ppt(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ppt_OpeningFcn, ...
                   'gui_OutputFcn',  @ppt_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

%executes when the gui opens
function ppt_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handlesArr = [handles.Flare_Angle_slider,handles.Flare_angle_edittext,handles.flare_x_edit,handles.flare_x_sld,handles.flare_angle_slider_large, handles.flare_angle_edit_large, handles.flare_length_slider_large, handles.flare_length_edit_large];
set(handlesArr, 'Enable', 'off');
handles1=[handles.micro_ppt,handles.large_ppt];
set(handles1,'BackgroundColor',[0.8,0.8,0.8]);
set(handles.main_menu,'BackgroundColor',[1,1,1]);
set(handles.Activate_Flared_Geometry_radiob,'value',0);
set(handles.flare_button_large,'value',0);
set(handles.main_panel,'Visible','On');
set(handles.micro_ppt_panel,'Visible','Off');
set(handles.large_ppt_panel,'Visible','Off');
guidata(hObject,handles);

function varargout = ppt_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function Discharge_energy_slider_Callback(hObject, eventdata, handles)
Den=get(hObject,'value');
set(handles.Discharge_energy_edittext,'string',num2str(Den));
guidata(hObject,handles);

function Discharge_energy_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Capacitance_slider_Callback(hObject, eventdata, handles)
Cap=get(hObject,'value');
set(handles.Capacitance_edittext,'string',num2str(Cap));
guidata(hObject,handles);

function Capacitance_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Pulse_frequency_slider_Callback(hObject, eventdata, handles)
Pfreq=get(hObject,'value');
set(handles.Pulse_frequency_edit,'string',num2str(Pfreq));
guidata(hObject,handles);

function Pulse_frequency_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%reset button for microppt
function Reset_push_Callback(hObject, eventdata, handles)
handles1=[handles.cap_inductance_sld,handles.Flare_Angle_slider,handles.prop_height_sld,handles.cap_inductance_sld,handles.Discharge_Energy_slider,handles.Pulse_frequency_slider,handles.electrode_Inductance_slider,handles.Capacitance_slider,handles.electrode_Resistance_slider,handles.ESR_slider];
set(handles1,'value',0);
handles2=[handles.Range_min_edit,handles.Range_max_edit,handles.prop_height_edit,handles.cap_inductance_edit,handles.Discharge_energy_edittext,handles.Pulse_frequency_edit,handles.electrode_inductance_edittext,handles.Capacitance_edittext,handles.electrode_Resistance_edit,handles.cap_inductance_edit,handles.ESR_edit];
set(handles2,'string','');
handles3=[handles.Propellant_width_edittext,handles.Flare_angle_edittext,handles.ed_spacing_edit,handles.ed_width_edit,handles.length_edit,handles.flare_x_edit,handles.prop_area_edit,handles.ar_edit,handles.EA_ratio_edittext];
set(handles3,'string','');
handles4=[handles.Resolution_edit,handles.time_results_edit,handles.thrust_edit,handles.Exhaust_velocity_edit,handles.Impulse_edit,handles.Specific_Impulse_edit,handles.Inductance_gradient_edit,handles.transfer_edit,handles.acceleration_edit,handles.Efficiency_edit,handles.Inductance_change_edit];
set(handles4,'string','');
handles5=[handles.Propellant_width_slider,handles.ed_spacing_slider,handles.ed_width_sld,handles.Electrode_length_slider,handles.flare_x_sld,handles.electrode_thickness_sld];
set(handles5,'value',0);
handles6=[handles.Electrode_shape_popup,handles.Dependent_Variable_popup,handles.Independent_Variable_popup,handles.Propellant_feeding_popup,handles.plot_choice];
set(handles6,'value',1);
handles7=[handles.spot_radius_micro,handles.spot_velocity,handles.ion_erosion,handles.Mean_ion_charge,handles.alpha_i,handles.current_per_spot,handles.electrode_thickness_edit,handles.dm_dt];
set(handles7,'string','');
axes(handles.axes1);
cla reset;
axes(handles.axes2);
cla reset;
axes(handles.axes3);
cla reset;
set(handles.Activate_Flared_Geometry_radiob,'value',0);

function Inductance_change_edit_Callback(hObject, eventdata, handles)

function Inductance_change_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function thrust_edit_Callback(hObject, eventdata, handles)

function thrust_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Exhaust_velocity_edit_Callback(hObject, eventdata, handles)

function Exhaust_velocity_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Impulse_edit_Callback(hObject, eventdata, handles)

function Impulse_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Efficiency_edit_Callback(hObject, eventdata, handles)

function Efficiency_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Inductance_gradient_edit_Callback(hObject, eventdata, handles)

function Inductance_gradient_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Specific_Impulse_edit_Callback(hObject, eventdata, handles)

function Specific_Impulse_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% plot function in micro ppt section
function Clear_plot2_push_Callback(hObject, eventdata, handles)
set(handles.Clear_plot2_push,'BackgroundColor',[1,1,1]);
axes(handles.axes2);
cla reset;
set(handles.Clear_plot2_push,'BackgroundColor',[0.85,0.85,0.85]);

function Save_plot2_push_Callback(hObject, eventdata, handles)
set(handles.Save_plot2_push,'BackgroundColor',[1,1,1]);
set(handles.axes2,'units','normalized');
fignew = figure('Visible','on'); % Invisible figure
newAxes = copyobj(handles.axes2,fignew); % Copy the appropriate axes
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
set(fignew,'CreateFcn','set(gcbf,''Visible'',''on'')'); % Make it visible upon loading
set(fignew, 'PaperPositionMode', 'auto');  
set(handles.Save_plot2_push,'BackgroundColor',[0.85,0.85,0.85]);

function Clear_plot3_push_Callback(hObject, eventdata, handles)
set(handles.Clear_plot3_push,'BackgroundColor',[1,1,1]);
axes(handles.axes3);
cla reset;
set(handles.Clear_plot3_push,'BackgroundColor',[0.85,0.85,0.85]);

function Save_plot3_push_Callback(hObject, eventdata, handles)
set(handles.Save_plot3_push,'BackgroundColor',[1,1,1]);
set(handles.axes3,'units','normalized');
fignew = figure('Visible','on'); % Invisible figure
newAxes = copyobj(handles.axes3,fignew); % Copy the appropriate axes
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
set(fignew,'CreateFcn','set(gcbf,''Visible'',''on'')'); % Make it visible upon loading
set(fignew, 'PaperPositionMode', 'auto');  
set(handles.Save_plot3_push,'BackgroundColor',[0.85,0.85,0.85]);

%plotting for microppt
function Plot_push_Callback(hObject, eventdata, handles)
set(handles.Plot_push,'BackgroundColor',[1,1,1]);
progress_f= waitbar(0,'Please wait...');
waitbar(.33,progress_f,'Loading your data');
pause(0.5);
switch get(handles.plot_choice,'Value')
    case 2
        axes(handles.axes2);
        cla reset;
    case 3
        axes(handles.axes3);
        cla reset;
end
min=(str2num(get(handles.Range_min_edit,'string')))*1e-3;
max=(str2num(get(handles.Range_max_edit,'string')))*1e-3;
res=(str2num(get(handles.Resolution_edit,'string')))*1e-3;
%calculations
digits(6)
u=1.2566370614*(1e-6);%permeability of free space
g=9.8;
E=str2num(get(handles.Discharge_energy_edittext,'string'));%discharge energy in J
f=str2num(get(handles.Pulse_frequency_edit,'string'));%pulse frequency in Hz
C=(str2num(get(handles.Capacitance_edittext,'string')))*1e-6;%capacitance in F
Le=(str2num(get(handles.electrode_inductance_edittext,'string')))*1e-9;%electrode inductance in H
Re=(str2num(get(handles.electrode_Resistance_edit,'string')))*1e-3;%electrode resistance in ohm                              
ESR=(str2num(get(handles.ESR_edit,'string')))*1e-3;%capacitor resistance in ohm 
Lcap=(str2num(get(handles.cap_inductance_edit,'string')))*1e-9;%capacitor inductance in H
h=(str2num(get(handles.ed_spacing_edit,'string')))*1e-3;%electrode spacing
w=(str2num(get(handles.ed_width_edit,'string')))*1e-3;%electrode width
l=(str2num(get(handles.length_edit,'string')))*1e-3;%electrode length
d=(str2num(get(handles.electrode_thickness_edit,'string')))*1e-3;%electrode thickness
t_span=(str2num(get(handles.time_results_edit,'string')))*1e-6;%time_span for plasma sheet to reach the end of electrodes                                
if (h/w >2)
    warndlg('It seems the aspect ratio is greater than 1. This makes your calculation into the catrgory of Large PPT. Kindly, check the dimensions again or use the Large PPT caculation page.');
else   
    L=Le+Lcap;%total initial resistance
    R=ESR+Re;%total initial resistance
    Vo=sqrt(2*E/C);%voltage
    zi=(R/2)*sqrt(C/L);%zeta
    wo=1/sqrt(L*C);
    syms t;%declaring t as a variaable
    if zi<1
        i2=@(t)((Vo^2.*exp(-2*zi*wo*t)*(sin(sqrt(1-zi^2)*wo*t))^2)/(L^2.*(1-zi^2).*wo^2));%current equation
        i=@(t)((Vo.*exp(-zi*wo*t)*(sin(sqrt(1-(zi^2)).*wo*t)))/(L.*sqrt(1-(zi^2)).*wo));
        v=@(t)((Vo.*zi.*exp(-zi*wo*t).*sin(sqrt(1-(zi^2)).*wo*t)./sqrt(1-(zi^2)))+(Vo.*exp(-zi*wo*t).*cos(sqrt(1-(zi^2)).*wo*t))); %voltage equation
    else
        if C==(4*L/(R^2))
            i2=@(t)(Vo^2/L^2)*(t^2.*exp(-2*wo*t));%current equation
            i=@(t)(Vo/L)*(t.*exp(-wo*t));
            v=@(t)(Vo.*exp(-wo*t).*(1+wo*t));%voltage equation
        else
            i2=@(t)(Vo^2.*(exp(-2*zi*wo*t))*(sinh(wo*t*sqrt(zi^2-1)))^2)/(L^2.*(zi^2-1).*wo^2);%current equation
            i=@(t)(Vo.*(exp(-zi*wo*t))*(sinh(wo*t*sqrt((zi^2)-1))))/(L.*sqrt((zi^2)-1).*wo);
            v=@(t)((Vo.*zi.*exp(-zi*wo*t).*sinh(sqrt((zi^2)-1).*wo*t)./sqrt((zi^2)-1))+(Vo.*exp(-zi*wo*t).*cosh(sqrt((zi^2)-1).*wo*t)));
        end
    end
    V=v(t_span);
    I=i(t_span);%loop current flowing through the circuit
    I_arc=Vo*sqrt(C/L);%arc current flowing through the circuit loop
    i2_integration=int(i2,t,[0,t_span]);%evaluation of integral of (i(t))^2
    Xp=vpa(i2_integration);%making the above answer definite
    %calculation of mbit
    me=9.1.*1e-31;%mass of electron
    e=1.6.*1e-19;%elementary charge
    eo=8.85.*1e-12;%Vacuum permittivity of free space(F/m)
    Ro=(1/3)*sqrt(C/L);
    Lp=0;%plasma inductance
    Te=11605;%temperature(in K)(1 eV)
    r=(str2num(get(handles.spot_radius_micro,'string')))*1e-3;%radius of the cathode spot
    S_spot=pi.*r^2;%Surface area of the cathode spot near the mixing region
    I_spot=(str2num(get(handles.current_per_spot,'string')));%For copper electrodes it has been experimentally observed that the current per observed spot
    So=S_spot*I_arc/I_spot;%Surface area of the initial plasma flow near the mixing region
    del_T= 92840;%temperature difference(8 eV)(in K)
    lambda=1.2*del_T;%latent heat of vaporization, specific heat=1.2
    alp_I=(str2num(get(handles.alpha_i,'string')));%Ion current normalised by arc current, the ratio between the ion current and the arc current in a cathodic plasma
    Q=(str2num(get(handles.Mean_ion_charge,'string')));%mean ion charge state number in plasma flow
    V_spot=(str2num(get(handles.spot_velocity,'string')))*1e3;%initial copper plasma has a velocity towards the direction of the anode
    Ne=I*alp_I./(Q*e*So*V_spot);%Electron density
    Ni=Ne/Q;%ion density
    A=23-log(Ne^0.5*Q*Te^-1.5);%Coulomb logarithm, factor by which small angle collisions are more effective than large angle collision
    neu=(3.62.*1e-6)*A*Ne*(Te^(-3/2))/Q;%The frequency at which electrons and ions collide
    sigma=(Ne.*(e^2))/(me*neu);%conductivity
    Rp=abs(V/I);%plasma resistance
    V_sheath=abs(V/2);
    d_sheath=sqrt(2*V_sheath*eo/(e*Ni));%thickness of the sheet
    waitbar(.67,progress_f,'Processing your data');
    pause(0.5);
    syms x;
    switch get(handles.Independent_Variable_popup,'Value')
        case 2
            h = min:res:max;
            %for non flared geometry
            if (get(handles.Activate_Flared_Geometry_radiob,'value'))==0
                m_bit=(sigma*d_sheath*d./(lambda.*h)).*(Rp^2+Lp^2/L*C)*(E/Ro);
                mbit=vpa(m_bit);
                switch get(handles.Electrode_shape_popup,'Value')
                    case 2
                        xp=(u.*h./(2.*mbit.*w)).*Xp;
                    case 3
                        tl=0.2*l;%part with electrode in tongue shape
                        xp1=(u.*h./(2.*mbit.*w)).*Xp;%calculation for rectangular part
                        xp1=vpa(xp1);
                        xp2=@(x)(((u.*h./(2.*mbit.*w)).*Xp)./(1-x/l));%calculation for tongue part
                        xp2=vpa(int(xp2,x,[0,tl]));
                        xp=xp1+xp2;
                end
            end
            %for flared geometry
            if (get(handles.Activate_Flared_Geometry_radiob,'value'))==1
                flare_l=((str2num(get(handles.flare_x_edit,'string')))*1e-3);
                alpha=str2num(get(handles.Flare_angle_edittext,'string'));%flare angle
                switch get(handles.Electrode_shape_popup,'Value')
                    case 2
                        m_bit1=(sigma*d_sheath*d./(lambda.*h)).*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for rectangular parallel part
                        m_bit1=vpa(m_bit1);
                        m_bit2=@(x)((sigma*d_sheath*d./(lambda.*(h+(tand(alpha))*x))).*(Rp^2+Lp^2./(L*C)).*(E/Ro));%calculation for rectangular flared part
                        m_bit2=vpa(int(m_bit2,x,[0,flare_l]));
                        mbit=m_bit1+m_bit2;
                        xp1=(u.*h./(2.*m_bit1.*w)).*Xp;%calculation for rectangular parallel part
                        xp1=vpa(xp1);
                        xp2=@(x)(((u.*(h+(tand(alpha))*x))./(2.*m_bit2*w)).*Xp);%calculation for rectangular flared part
                        xp2=vpa(int(xp2,x,[0,flare_l]));
                        xp=xp1+xp2;
                    case 3
                        m_bit1=(sigma*d_sheath*d./(lambda.*h))*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for rectangular parallel part
                        m_bit1=vpa(m_bit1);
                        m_bit2=@(x)((sigma*d_sheath.*d)./(lambda.*(h+(tand(alpha))*x))).*(Rp^2+Lp^2/L*C)*(E/Ro);  %calculation for flared part
                        m_bit2=vpa(int(m_bit2,x,[0,flare_l]));
                        mbit=m_bit1+m_bit2;
                        xp1=(u.*h./(2.*m_bit1.*w)).*Xp;%calculation for rectangular parallel part
                        flare_l_a=0.8.*flare_l;%part with electrode in rectangular shape
                        flare_l_b=0.2.*flare_l;%part with electrode in tongue shape
                        xp2=@(x)(((u.*(h+(tand(alpha))*x))./(2.*m_bit2.*w)).*Xp);
                        xp2=vpa(int(xp2,x,[0,flare_l_a]));
                        xp3=@(x)(((u.*(h+(tand(alpha))*x))./(2.*m_bit2.*(w*(1-x/l)))).*Xp);
                        xp3=vpa(int(xp3,x,[0,flare_l_b]));
                        xp=xp1+xp2+xp3;
                end
            end
            ve=xp;%exhaust velocity
            Isp = ve./g;%Specific Impulse
            impulse_bit=Isp.*mbit.*g*f; %impulse bit
            thrust=impulse_bit.*f;  %thrust 
            Efficiency=mbit.*(ve.^2)./(2*E);%efficiency
            h=h*1e3;
            switch get(handles.Dependent_Variable_popup,'Value')
                case 2
                    plot(h,Efficiency*1e2,'-*');
                    xlabel('electrode spacing(mm)');
                    ylabel('Thruster Efficiency(%)');
                case 3
                    plot(h,impulse_bit*1e6,'-*');
                    xlabel('electrode spacing(mm)');
                    ylabel('Impulse Bit(uNs)');
                case 4
                    plot(h,ve,'-*');
                    xlabel('electrode spacing(mm)');
                    ylabel('Exhaust Velocity(m/s)');
                case 5
                    plot(h,thrust*1e6,'-*');
                    xlabel('electrode spacing(mm)');
                    ylabel('Thrust(uN)');
            end
        case 3
            w = min:res:max;
            %for non flared geometry
            if (get(handles.Activate_Flared_Geometry_radiob,'value'))==0
                m_bit=(sigma*d_sheath*d./(lambda.*h)).*(Rp^2+Lp^2/L*C)*(E/Ro);
                mbit=vpa(m_bit);
                switch get(handles.Electrode_shape_popup,'Value')
                    case 2
                        xp=(u.*h./(2.*mbit.*w)).*Xp;
                    case 3
                        tl=0.2*l;%part with electrode in tongue shape
                        xp1=(u.*h./(2.*mbit.*w)).*Xp;%calculation for rectangular part
                        xp1=vpa(xp1);
                        xp2=@(x)(((u.*h./(2.*mbit.*w)).*Xp)./(1-x/l));%calculation for tongue part
                        xp2=vpa(int(xp2,x,[0,tl]));
                        xp=xp1+xp2;
                end
            end
            %for flared geometry
            if (get(handles.Activate_Flared_Geometry_radiob,'value'))==1
                flare_l=((str2num(get(handles.flare_x_edit,'string')))*1e-3);
                alpha=str2num(get(handles.Flare_angle_edittext,'string'));%flare angle
                switch get(handles.Electrode_shape_popup,'Value')
                    case 2
                        m_bit1=(sigma*d_sheath*d./(lambda.*h)).*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for rectangular parallel part
                        m_bit1=vpa(m_bit1);
                        m_bit2=@(x)((sigma*d_sheath*d./(lambda.*(h+(tand(alpha))*x))).*(Rp^2+Lp^2./(L*C)).*(E/Ro));%calculation for rectangular flared part
                        m_bit2=vpa(int(m_bit2,x,[0,flare_l]));
                        mbit=m_bit1+m_bit2;
                        xp1=(u.*h./(2.*m_bit1.*w)).*Xp;%calculation for rectangular parallel part
                        xp1=vpa(xp1);
                        xp2=@(x)(((u.*(h+(tand(alpha))*x))./(2.*m_bit2*w)).*Xp);%calculation for rectangular flared part
                        xp2=vpa(int(xp2,x,[0,flare_l]));
                        xp=xp1+xp2;
                    case 3
                        m_bit1=(sigma*d_sheath*d./(lambda.*h))*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for rectangular parallel part
                        m_bit1=vpa(m_bit1);
                        m_bit2=@(x)((sigma*d_sheath.*d)./(lambda.*(h+(tand(alpha))*x))).*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for flared part
                        m_bit2=vpa(int(m_bit2,x,[0,flare_l]));
                        mbit=m_bit1+m_bit2;
                        xp1=(u.*h./(2.*m_bit1.*w)).*Xp;%calculation for rectangular parallel part
                        flare_l_a=0.8.*flare_l;%part with electrode in rectangular shape
                        flare_l_b=0.2.*flare_l;%part with electrode in tongue shape
                        xp2=@(x)(((u.*(h+(tand(alpha))*x))./(2.*m_bit2.*w)).*Xp);
                        xp2=vpa(int(xp2,x,[0,flare_l_a]));
                        xp3=@(x)(((u.*(h+(tand(alpha))*x))./(2.*m_bit2.*(w*(1-x/l)))).*Xp);
                        xp3=vpa(int(xp3,x,[0,flare_l_b]));
                        xp=xp1+xp2+xp3;
                end  
            end
            ve=xp;%exhaust velocity
            Isp = ve./g;%Specific Impulse
            impulse_bit=Isp*mbit*g*f;%impulse bit
            thrust=impulse_bit*f;%thrust 
            Efficiency=mbit.*(ve.^2)./(2*E);%efficiency
            w=w.*1e3;
            switch get(handles.Dependent_Variable_popup,'Value')
                case 2
                    plot(w,Efficiency*1e2,'-*');
                    xlabel('electrode width(mm)');
                    ylabel('Thruster Efficiency(%)');
                case 3
                    plot(w,impulse_bit*1e6,'-*');
                    xlabel('electrode width(mm)');
                    ylabel('Impulse Bit(uNs)');
                case 4
                    plot(w,ve,'-*');
                    xlabel('electrode width(mm)');
                    ylabel('Exhaust Velocity(m/s)');
                case 5
                    plot(w,thrust*1e6,'-*');
                    xlabel('electrode width(mm)');
                    ylabel('Thrust(uN)');
            end
    end
end
waitbar(1,progress_f,'Finished');
grid on;
pause(0.5);
close(progress_f);
set(handles.Plot_push,'BackgroundColor',[0.85,0.85,0.85]);
guidata(hObject,handles);

function Dependent_Variable_popup_Callback(hObject, eventdata, handles)

function Dependent_Variable_popup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Independent_Variable_popup_Callback(hObject, eventdata, handles)

function Independent_Variable_popup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Range_min_edit_Callback(hObject, eventdata, handles)

function Range_min_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Range_max_edit_Callback(hObject, eventdata, handles)

function Range_max_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Resolution_edit_Callback(hObject, eventdata, handles)

function Resolution_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Propellant_feeding_popup_Callback(hObject, eventdata, handles)
   
function Propellant_feeding_popup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Electrode_shape_popup_Callback(hObject, eventdata, handles)

function Electrode_shape_popup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Propellant_width_slider_Callback(hObject, eventdata, handles)
width=get(hObject,'value');
set(handles.Propellant_width_edittext,'string',num2str(width));
guidata(hObject,handles);

function Propellant_width_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Activate_Flared_Geometry_radiob_Callback(hObject, eventdata, handles)
handlesArray = [handles.Flare_Angle_slider,handles.Flare_angle_edittext,handles.flare_x_sld,handles.flare_x_edit];
if get(hObject,'Value') == 1
    set(handlesArray, 'Enable', 'on');
else
    set(handlesArray, 'Enable', 'off');
end

function Flare_Angle_slider_Callback(hObject, eventdata, handles)
fangle=get(hObject,'value');
set(handles.Flare_angle_edittext,'string',num2str(fangle));
guidata(hObject,handles);

function Flare_Angle_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Discharge_energy_edittext_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.Discharge_Energy_slider,'value',str2num(edit));
energy = str2num(get(handles.Discharge_energy_edittext,'String'));
minSliderValue = get(handles.Discharge_Energy_slider, 'Min');
maxSliderValue = get(handles.Discharge_Energy_slider, 'Max');
if energy < minSliderValue
    e = minSliderValue ;
    set(handles.Discharge_energy_edittext,'String',num2str(e));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.Discharge_Energy_slider, 'Value', e);
elseif  energy > maxSliderValue 
    e = maxSliderValue;
    set(handles.Discharge_energy_edittext,'String',num2str(e));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.Discharge_Energy_slider, 'Value', e);
end
guidata(hObject,handles);

function Discharge_energy_edittext_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_inductance_edittext_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.electrode_Inductance_slider,'value',str2num(edit));
ind = str2num(get(handles.electrode_inductance_edittext,'String'));
minSliderValue = get(handles.electrode_Inductance_slider, 'Min');
maxSliderValue = get(handles.electrode_Inductance_slider, 'Max');
if ind < minSliderValue
    in = minSliderValue;
    set(handles.electrode_inductance_edittext,'String',num2str(in));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_Inductance_slider, 'Value', in);
elseif  ind > maxSliderValue 
    in = maxSliderValue;
    set(handles.electrode_inductance_edittext,'String',num2str(in));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_Inductance_slider, 'Value', in);
end
guidata(hObject,handles);

function electrode_inductance_edittext_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Capacitance_edittext_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.Capacitance_slider,'value',str2num(edit));
cap = str2num(get(handles.Capacitance_edittext,'String'));
minSliderValue = get(handles.Capacitance_slider, 'Min');
maxSliderValue = get(handles.Capacitance_slider, 'Max');
if cap < minSliderValue
    c = minSliderValue ;
    set(handles.Capacitance_edittext,'String',num2str(c));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.Capacitance_slider, 'Value', c);
elseif  cap > maxSliderValue 
    c = maxSliderValue;
    set(handles.Capacitance_edittext,'String',num2str(c));
    warndlg('Value Out of Range! Please Enter Again');
        set(handles.Capacitance_slider, 'Value', c);
end
guidata(hObject,handles);

function Capacitance_edittext_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Propellant_width_edittext_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.Propellant_width_slider,'value',str2num(edit));
width = str2num(get(handles.Propellant_width_edittext,'String'));
minSliderValue = get(handles.Propellant_width_slider, 'Min');
maxSliderValue = get(handles.Propellant_width_slider, 'Max');
if width < minSliderValue
    w = minSliderValue ;
    set(handles.Propellant_width_edittext,'String',num2str(w));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.Propellant_width_slider, 'Value', w);
elseif  width > maxSliderValue 
    w = maxSliderValue;
    set(handles.Propellant_width_edittext,'String',num2str(w));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.Propellant_width_slider, 'Value', w);
end
guidata(hObject,handles);

function Propellant_width_edittext_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%final results calculation for micro ppt
function Calculate_results_push_Callback(hObject, eventdata, handles)
set(handles.Calculate_results_push,'BackgroundColor',[1,1,1]);
digits(7);
progress_f= waitbar(0,'Please wait...');
waitbar(.33,progress_f,'Loading your data');
pause(0.5);
%plot Guman and Paccani Graph
axes(handles.axes1);
cla reset;
YISP= [90 200 290 400 460 600 710 805 900 990 1070 1140 1215];
XEA = [0.2 0.6 1 1.6 2 3 4 5 6 7 8 9 10];
X_E_A = linspace(min(XEA), max(XEA), 150);% Evenly-Spaced Interpolation Vector for guman
Y_ISP = interp1(XEA, YISP, X_E_A, 'spline', 'extrap');
p1=plot(X_E_A, Y_ISP,'LineWidth',0.05);
xlabel('E/A (J/cm^2)');
ylabel('Predicted Isp (s)');
grid on;
hold on;
YISP1= [200 260 340 380 510 600 666 725 777 835 880 955 1025 1090 1150 1210 1260];
XEA1 = [0.1 0.2 0.3 0.4 1 1.5 2 2.5 3 3.5 4 5 6 7 8 9 10];
X_E_A1 = linspace(min(XEA1), max(XEA1), 150);% Evenly-Spaced Interpolation Vector for peccani
Y_ISP1 = interp1(XEA1, YISP1, X_E_A1, 'spline', 'extrap');
p2=plot(X_E_A1, Y_ISP1,'LineWidth',0.05);
legend([p1(1),p2(1)],'Guman','Gessini-Paccani','Location','Northwest');
%calculations
u=1.2566370614*(1e-6);%permeability of free space
g=9.8;
E=str2num(get(handles.Discharge_energy_edittext,'string'));%discharge energy in J
f=str2num(get(handles.Pulse_frequency_edit,'string'));%pulse frequency in Hz
C=(str2num(get(handles.Capacitance_edittext,'string')))*1e-6;%capacitance in F
Le=(str2num(get(handles.electrode_inductance_edittext,'string')))*1e-9;%electrode inductance in H
Re=(str2num(get(handles.electrode_Resistance_edit,'string')))*1e-3; %electrode resistance in ohm                              
ESR=(str2num(get(handles.ESR_edit,'string')))*1e-3;%capacitor resistance in ohm 
Lcap=(str2num(get(handles.cap_inductance_edit,'string')))*1e-9; %capacitor inductance in H
h=(str2num(get(handles.ed_spacing_edit,'string')))*1e-3;%electrode spacing
w=(str2num(get(handles.ed_width_edit,'string')))*1e-3; %electrode width
l=(str2num(get(handles.length_edit,'string')))*1e-3; %electrode length
d=(str2num(get(handles.electrode_thickness_edit,'string')))*1e-3; %electrode thickness
t_span=(str2num(get(handles.time_results_edit,'string')))*1e-6;%time_span for plasma sheet to reach the end of electrodes                                
if (h/w >2)
    warndlg('It seems the aspect ratio is greater than 1. This makes your calculation into the catrgory of Large PPT. Kindly, check the dimensions again or use the Large PPT caculation page.');
else  
    L=Le+Lcap;%total initial resistance
    R=ESR+Re;%total initial resistance
    Vo=sqrt(2*E/C);%voltage
    zi=(R/2)*sqrt(C/L);%zeta
    wo=1/sqrt(L*C);
    syms t;%declaring t as a variaable
    if zi<1
        i2=@(t)((Vo^2.*exp(-2*zi*wo*t)*(sin(sqrt(1-zi^2)*wo*t))^2)/(L^2.*(1-zi^2).*wo^2)); %current equation
        i=@(t)((Vo.*exp(-zi*wo*t)*(sin(sqrt(1-(zi^2)).*wo*t)))/(L.*sqrt(1-(zi^2)).*wo));
        v=@(t)((Vo.*zi.*exp(-zi*wo*t).*sin(sqrt(1-(zi^2)).*wo*t)./sqrt(1-(zi^2)))+(Vo.*exp(-zi*wo*t).*cos(sqrt(1-(zi^2)).*wo*t))); %voltage equation
    else
        if C==(4*L/(R^2))
            i2=@(t)(Vo^2/L^2)*(t^2.*exp(-2*wo*t));%current equation
            i=@(t)(Vo/L)*(t.*exp(-wo*t));
            v=@(t)(Vo.*exp(-wo*t).*(1+wo*t));%voltage equation
        else
            i2=@(t)(Vo^2.*(exp(-2*zi*wo*t))*(sinh(wo*t*sqrt(zi^2-1)))^2)/(L^2.*(zi^2-1).*wo^2);%current equation
            i=@(t)(Vo.*(exp(-zi*wo*t))*(sinh(wo*t*sqrt((zi^2)-1))))/(L.*sqrt((zi^2)-1).*wo);
            v=@(t)((Vo.*zi.*exp(-zi*wo*t).*sinh(sqrt((zi^2)-1).*wo*t)./sqrt((zi^2)-1))+(Vo.*exp(-zi*wo*t).*cosh(sqrt((zi^2)-1).*wo*t)));
        end
    end
    V=v(t_span);
    I=i(t_span); %loop current flowing through the circuit
    I_arc=Vo*sqrt(C/L);%arc current flowing through the circuit loop
    i2_integration=int(i2,t,[0,t_span]);%evaluation of integral of (i(t))^2
    Xp=vpa(i2_integration);%making the above answer definite
    %calculation of mbit
    me=9.1.*1e-31;%mass of electron
    e=1.6.*1e-19;%elementary charge
    eo=8.85.*1e-12;%Vacuum permittivity of free space(F/m)
    Ro=(1/3)*sqrt(C/L);
    Lp=0; %plasma inductance
    Te=11605; %temperature(in K)(1 eV)
    r=(str2num(get(handles.spot_radius_micro,'string')))*1e-3;%radius of the cathode spot
    S_spot=pi.*r^2;%Surface area of the cathode spot near the mixing region
    I_spot=(str2num(get(handles.current_per_spot,'string')));%For copper electrodes it has been experimentally observed that the current per observed spot
    So=S_spot*I_arc/I_spot;%Surface area of the initial plasma flow near the mixing region
    del_T= 92840;%temperature difference(8 eV)(in K)
    lambda=1.2*del_T;%latent heat of vaporization, specific heat=1.2
    alp_I=(str2num(get(handles.alpha_i,'string')));%Ion current normalised by arc current, the ratio between the ion current and the arc current in a cathodic plasma
    Q=(str2num(get(handles.Mean_ion_charge,'string')));%mean ion charge state number in plasma flow
    V_spot=(str2num(get(handles.spot_velocity,'string')))*1e3;%initial copper plasma has a velocity towards the direction of the anode
    Ne=I*alp_I./(Q*e*So*V_spot);%Electron density
    Ni=Ne/Q;%ion density
    A=23-log(Ne^0.5*Q*Te^-1.5);%Coulomb logarithm, factor by which small angle collisions are more effective than large angle collision
    neu=(3.62.*1e-6)*A*Ne*(Te^(-3/2))/Q;%The frequency at which electrons and ions collide
    sigma=(Ne.*(e^2))/(me*neu);%conductivity
    Rp=abs(V/I);%plasma resistance
    V_sheath=abs(V/2);
    d_sheath=sqrt(2*V_sheath*eo/(e*Ni));%thickness of the sheet
    I_sheath=pi*(Re^2)*(2.33e-6).*((V_sheath^1.5)/(d_sheath^2));
    R_sheath=V_sheath/I_sheath;%sheath resistance
    Ztot=R+Rp+R_sheath;%total resistance of the circuit
    waitbar(.67,progress_f,'Processing your data');
    pause(0.5);
    syms x;
    %for non flared geometry
    if (get(handles.Activate_Flared_Geometry_radiob,'value'))==0
        m_bit=(sigma*d_sheath*d/(lambda*h))*(Rp^2+Lp^2/L*C)*(E/Ro);%m_ablated = m_bit
        mbit=vpa(m_bit);
        switch get(handles.Electrode_shape_popup,'Value')
            case 2
                xp=(u*h/(2.*mbit.*w)).*Xp;
                ind_grad=u.*h/w;
            case 3
                tl=0.2*l;%part with electrode in tongue shape
                xp1=(u*h/(2.*mbit.*w)).*Xp;%calculation for rectangular part
                xp1=vpa(xp1);
                xp2=@(x)((u*h/(2*mbit*w.*(1-x/l)))*Xp);%calculation for tongue part
                xp2=vpa(int(xp2,x,[0,tl]));
                xp=xp1+xp2;
                ind_grad=u.*h/w;
        end
        ve=xp; %exhaust velocity
        Isp = ve/g; %Specific Impulse
        Efficiency=mbit.*ve^2/(2*E);%efficiency
        impulse_bit=Isp*mbit*g*f; %impulse bit
        thrust=impulse_bit*f;  %thrust
        tr_eff=1-(ESR/Ztot);%transfer efficiency
        acc_eff=Efficiency/tr_eff;%acceleration efficiency
        hprop=(str2num(get(handles.prop_height_edit,'string')))*1e-3;
        wprop=(str2num(get(handles.Propellant_width_edittext,'string')))*1e-3;
        aspect_ratio=h/w;
        set(handles.ar_edit,'string',num2str(aspect_ratio));
        switch get(handles.Propellant_feeding_popup,'Value')
            case 2
                n=1;
            case 3
                n=2;
            case 4
                n=3;
        end
        expo_area=(n*hprop*wprop)*(1e4);
        E_A=E./expo_area;
        gamma=(str2num(get(handles.ion_erosion,'string')));%ion erosion rate(in µg/C)
        %considering the number of pulses to be 5e5
        dm_dt=I.*gamma/5e5;%rate of mass loss from the electrode surface as a function of the discharge current
        set(handles.dm_dt,'string',num2str(dm_dt));
        set(handles.prop_area_edit,'string',num2str(expo_area));
        set(handles.EA_ratio_edittext,'string',num2str(E_A));
        set(handles.Specific_Impulse_edit,'string',char(Isp));
        set(handles.Exhaust_velocity_edit,'string',char(ve));
        set(handles.Inductance_gradient_edit,'string',num2str(ind_grad*1e6));
        set(handles.Impulse_edit,'string',char(impulse_bit*1e6));
        set(handles.Efficiency_edit,'string',char(Efficiency*1e2));
        set(handles.thrust_edit,'string',char(thrust*1e6));
        set(handles.transfer_edit,'string',num2str(tr_eff*1e2));
        set(handles.acceleration_edit,'string',char(acc_eff*1e2));
        set(handles.Inductance_change_edit,'string','nan');
    end
    %for flared geometry
    if (get(handles.Activate_Flared_Geometry_radiob,'value'))==1
        flare_l= ((str2num(get(handles.flare_x_edit,'string')))*1e-3);
        alpha=str2num(get(handles.Flare_angle_edittext,'string'));%flare angle
        switch get(handles.Electrode_shape_popup,'Value')
            case 2
                m_bit1=(sigma*d_sheath*d/(lambda*h))*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for rectangular parallel part
                m_bit1=vpa(m_bit1);
                m_bit2=@(x)((sigma*d_sheath*d/(lambda.*(h+(tand(alpha))*x))).*(Rp^2+Lp^2./(L*C)).*(E/Ro));%calculation for rectangular flared part
                m_bit2=vpa(int(m_bit2,x,[0,flare_l]));
                mbit=m_bit1+m_bit2;
                xp1=(u*h/(2*m_bit1*w))*Xp;%calculation for rectangular parallel part
                xp1=vpa(xp1);
                xp2=@(x)(((u.*(h+(tand(alpha))*x))/(2*m_bit2*w))*Xp);%calculation for rectangular flared part
                xp2=vpa(int(xp2,x,[0,flare_l]));
                xp=xp1+xp2;
                ind_grad1=u.*h/w;%calculation for rectangular parallel part
                ind_grad2=@(x)(u.*((h+tand(alpha)*x))/w);
                ind_grad2=vpa(int(ind_grad2,x,[0,flare_l]));%calculation for rectangular flared part
                ind_grad=ind_grad1+ind_grad2;
                %calculation of inductance change
                fun_x=@(x)(u/(2*pi))*((2*(h+(tand(alpha))*x)/w)*(pi-2*atan((h+tand(alpha)*x)/w)+(h+tand(alpha)*x)/w.*log((h+tand(alpha)*x)))-2*log(w)+(1-((h+tand(alpha)*x)/w).^2)*(log((h+tand(alpha)*x).^2+w.^2)));
            case 3
                m_bit1=(sigma*d_sheath*d/(lambda*h))*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for rectangular parallel part
                m_bit1=vpa(m_bit1);
                m_bit2=@(x)((sigma*d_sheath.*d)./(lambda.*(h+(tand(alpha))*x)))*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for flared part
                m_bit2=vpa(int(m_bit2,x,[0,flare_l]));
                mbit=m_bit1+m_bit2;
                xp1=(u*h/(2*m_bit1.*w))*Xp;%calculation for rectangular parallel part
                flare_l_a=0.8.*flare_l;%part with electrode in rectangular shape
                flare_l_b=0.2.*flare_l;%part with electrode in tongue shape
                xp2=@(x)(((u.*(h+(tand(alpha))*x))/(2*m_bit2.*w))*Xp);
                xp2=vpa(int(xp2,x,[0,flare_l_a]));
                xp3=@(x)(((u.*(h+(tand(alpha))*x))/(2*m_bit2.*(w*(1-x/l))))*Xp);
                xp3=vpa(int(xp3,x,[0,flare_l_b]));
                xp=xp1+xp2+xp3;
                ind_grad1=u.*h/w;%calculation for rectangular parallel part
                ind_grad2=@(x)(u*((h+tand(alpha)*x))/w);%calculation for rectangular flared part
                ind_grad2=vpa(int(ind_grad2,x,[0,flare_l_a]));
                ind_grad3=@(x)(u*((h+tand(alpha)*x))/(w*(1-x/l)));%calculation for tongue flared part
                ind_grad3=vpa(int(ind_grad3,x,[0,flare_l_b]));
                ind_grad=ind_grad1+ind_grad2+ind_grad3;
                %calculation of inductance change
                fun_x=@(x)((u/(2*pi))*((2*(h+tand(alpha)*x)/(w*(1-x/l)))*(pi-2*atan((h+tand(alpha)*x)/(w*(1-x/l)))+(h+tand(alpha)*x)/(w*(1-x/l)).*log((h+tand(alpha)*x)))-2*log((w*(1-x/l)))+(1-((h+tand(alpha)*x)/(w*(1-x/l))).^2)*(log((h+tand(alpha)*x).^2+(w*(1-x/l)).^2))));
        end
        del_L=int(fun_x,x,[0,l]);
        delL=vpa(del_L,5);
        ve=xp;%exhaust velocity
        Isp = ve/g;%Specific Impulse
        Efficiency=mbit.*ve^2/(2*E);%efficiency
        impulse_bit=Isp*mbit*g*f;
        thrust=impulse_bit*f;
        tr_eff=1-(ESR/Ztot);%transfer efficiency
        acc_eff=Efficiency/tr_eff ;%acceleration efficiency
        hprop=(str2num(get(handles.prop_height_edit,'string')))*1e-3;
        wprop=(str2num(get(handles.Propellant_width_edittext,'string')))*1e-3;
        aspect_ratio=h/w;
        set(handles.ar_edit,'string',num2str(aspect_ratio));
        switch get(handles.Propellant_feeding_popup,'Value')
            case 2
                n=1;
            case 3
                n=2;
            case 4
                n=3;
        end
        expo_area=(n*hprop*wprop)*(1e4);
        E_A=E./expo_area;
        gamma=(str2num(get(handles.ion_erosion,'string')));%ion erosion rate
        %considering the number of pulses to be 5e5
        dm_dt=I.*gamma/5e5;%rate of mass loss from the electrode surface as a function of the discharge current
        set(handles.dm_dt,'string',num2str(dm_dt));
        set(handles.prop_area_edit,'string',num2str(expo_area));
        set(handles.EA_ratio_edittext,'string',num2str(E_A));
        set(handles.Specific_Impulse_edit,'string',char(Isp));
        set(handles.Exhaust_velocity_edit,'string',char(ve));
        set(handles.Inductance_gradient_edit,'string',char(ind_grad*1e6));
        set(handles.Impulse_edit,'string',char(impulse_bit*1e6));
        set(handles.Efficiency_edit,'string',char(Efficiency*1e2));
        set(handles.thrust_edit,'string',char(thrust*1e6));
        set(handles.transfer_edit,'string',num2str(tr_eff*1e2));
        set(handles.acceleration_edit,'string',char(acc_eff*1e2));
        set(handles.Inductance_change_edit,'string',char(delL*1e9));
    end
end
waitbar(1,progress_f,'Finished');
pause(0.5);
close(progress_f);
set(handles.Calculate_results_push,'BackgroundColor',[0.85,0.85,0.85]);
guidata(hObject,handles);

function time_results_edit_Callback(hObject, eventdata, handles)

function time_results_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Flare_angle_edittext_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.Flare_Angle_slider,'value',str2num(edit));
angle = str2num(get(handles.Flare_angle_edittext,'String'));
minSliderValue = get(handles.Flare_Angle_slider, 'Min');
maxSliderValue = get(handles.Flare_Angle_slider, 'Max');
if angle < minSliderValue
    a = minSliderValue;
    set(handles.Flare_angle_edittext,'String',num2str(a));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.Flare_Angle_slider, 'Value', a);
elseif  angle > maxSliderValue 
    a = maxSliderValue;
    set(handles.Flare_angle_edittext,'String',num2str(a));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.Flare_Angle_slider, 'Value', a);
end
guidata(hObject,handles);

function Flare_angle_edittext_CreateFcn(hObject, eventdata, handles) 
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Discharge_Energy_slider_Callback(hObject, eventdata, handles)
Den=get(hObject,'value');
set(handles.Discharge_energy_edittext,'string',num2str(Den));
guidata(hObject,handles);

function Discharge_Energy_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function electrode_Inductance_slider_Callback(hObject, eventdata, handles)
Inductance=get(hObject,'value');
set(handles.electrode_inductance_edittext,'string',num2str(Inductance));
guidata(hObject,handles);

function electrode_Inductance_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function electrode_Resistance_edit_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.electrode_Resistance_slider,'value',str2num(edit));
res = str2num(get(handles.electrode_Resistance_edit,'String'));
minSliderValue = get(handles.electrode_Resistance_slider, 'Min');
maxSliderValue = get(handles.electrode_Resistance_slider, 'Max');
if res < minSliderValue
     r = minSliderValue ;
    set(handles.electrode_Resistance_edit,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
   set(handles.electrode_Resistance_slider, 'Value', r);
elseif  res > maxSliderValue 
    r = maxSliderValue;
    set(handles.electrode_Resistance_edit,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_Resistance_slider, 'Value', r);
end
guidata(hObject,handles);

function electrode_Resistance_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Pulse_frequency_edit_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.Pulse_frequency_slider,'value',str2num(edit));
freq = str2num(get(handles.Pulse_frequency_edit,'String'));
minSliderValue = get(handles.Pulse_frequency_slider, 'Min');
maxSliderValue = get(handles.Pulse_frequency_slider, 'Max');
if freq < minSliderValue
     f = minSliderValue ;
    set(handles.Pulse_frequency_edit,'String',num2str(f));
    warndlg('Value Out of Range! Please Enter Again');
   set(handles.Pulse_frequency_slider, 'Value', f);
elseif  freq > maxSliderValue 
    f = maxSliderValue;
    set(handles.Pulse_frequency_edit,'String',num2str(f));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.Pulse_frequency_slider, 'Value', f);
end
guidata(hObject,handles);

function Pulse_frequency_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_Resistance_slider_Callback(hObject, eventdata, handles)
resistance=get(hObject,'value');
set(handles.electrode_Resistance_edit,'string',num2str(resistance));
guidata(hObject,handles);

function electrode_Resistance_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function EA_ratio_edittext_Callback(hObject, eventdata, handles)

function EA_ratio_edittext_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ar_edit_Callback(hObject, eventdata, handles)

function ar_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function prop_area_edit_Callback(hObject, eventdata, handles)

function prop_area_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    
function ESR_slider_Callback(hObject, eventdata, handles)
esr=get(hObject,'value');
set(handles.ESR_edit,'string',num2str(esr));
guidata(hObject,handles);

function ESR_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function ESR_edit_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.ESR_slider,'value',str2num(edit));
res = str2num(get(handles.ESR_edit,'String'));
minSliderValue = get(handles.ESR_slider, 'Min');
maxSliderValue = get(handles.ESR_slider, 'Max');
if res < minSliderValue
    r = minSliderValue ;
    set(handles.ESR_edit,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.ESR_slider, 'Value', r);
elseif  res > maxSliderValue
     r = maxSliderValue;
    set(handles.ESR_edit,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
   set(handles.ESR_slider, 'Value', r);
end
guidata(hObject,handles);

function ESR_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function transfer_edit_Callback(hObject, eventdata, handles)

function transfer_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function acceleration_edit_Callback(hObject, eventdata, handles)

function acceleration_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ztot_slider_Callback(hObject, eventdata, handles)
ztot=get(hObject,'value');
set(handles.Ztot_edit,'string',num2str(ztot));
guidata(hObject,handles);

function Ztot_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Ztot_edit_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.Ztot_slider,'value',str2num(edit));
res = str2num(get(handles.Ztot_edit,'String'));
minSliderValue = get(handles.Ztot_slider, 'Min');
maxSliderValue = get(handles.Ztot_slider, 'Max');
if res < minSliderValue
    r = minSliderValue;
    set(handles.Ztot_edit,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.Ztot_slider, 'Value', r);
elseif  res > maxSliderValue 
    r = maxSliderValue;
    set(handles.Ztot_edit,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.Ztot_slider, 'Value', r);
end
guidata(hObject,handles);

function Ztot_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_spacing_edit_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.ed_spacing_slider,'value',str2num(edit));
x = str2num(get(handles.ed_spacing_edit,'String'));
minSliderValue = get(handles.ed_spacing_slider, 'Min');
maxSliderValue = get(handles.ed_spacing_slider, 'Max');
if x < minSliderValue
    s = minSliderValue ;
    set(handles.ed_spacing_edit,'String',num2str(s));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.ed_spacing_slider, 'Value', s);
    elseif  x > maxSliderValue
    s = maxSliderValue;
    set(handles.ed_spacing_edit,'String',num2str(s));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.ed_spacing_slider, 'Value', s);
end
guidata(hObject,handles);

function ed_spacing_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_width_edit_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.ed_width_sld,'value',str2num(edit));
x = str2num(get(handles.ed_width_edit,'String'));
minSliderValue = get(handles.ed_width_sld, 'Min');
maxSliderValue = get(handles.ed_width_sld, 'Max');
if x < minSliderValue
     r = minSliderValue;
    set(handles.ed_width_edit,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
   set(handles.ed_width_sld, 'Value', r);
elseif  x > maxSliderValue 
    r = maxSliderValue;
    set(handles.ed_width_edit,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.ed_width_sld, 'Value', r);
end
guidata(hObject,handles);

function ed_width_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_spacing_slider_Callback(hObject, eventdata, handles)
spacing=get(hObject,'value');
set(handles.ed_spacing_edit,'string',num2str(spacing));
guidata(hObject,handles);

function ed_spacing_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function flare_x_edit_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.flare_x_sld,'value',str2num(edit));
x = str2num(get(handles.flare_x_edit,'String'));
minSliderValue = get(handles.flare_x_sld, 'Min');
maxSliderValue = get(handles.flare_x_sld, 'Max');
if x < minSliderValue
    r = minSliderValue ;
    set(handles.flare_x_edit,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.flare_x_sld, 'Value', r);
elseif  x > maxSliderValue 
     r = maxSliderValue;
    set(handles.flare_x_edit,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
   set(handles.flare_x_sld, 'Value', r);
end
guidata(hObject,handles);

function flare_x_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function flare_x_sld_Callback(hObject, eventdata, handles)
width=get(hObject,'value');
set(handles.flare_x_edit,'string',num2str(width));
guidata(hObject,handles);

function flare_x_sld_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function ed_width_sld_Callback(hObject, eventdata, handles)
width=get(hObject,'value');
set(handles.ed_width_edit,'string',num2str(width));
guidata(hObject,handles);

function ed_width_sld_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function length_edit_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.Electrode_length_slider,'value',str2num(edit));
len = str2num(get(handles.length_edit,'String'));
minSliderValue = get(handles.Electrode_length_slider, 'Min');
maxSliderValue = get(handles.Electrode_length_slider, 'Max');
if len < minSliderValue
    l = minSliderValue;
    set(handles.length_edit,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.Electrode_length_slider, 'Value', l);
elseif  len > maxSliderValue 
     l = maxSliderValue;
    set(handles.length_edit,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.Electrode_length_slider, 'Value', l);
end
guidata(hObject,handles);

function length_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Electrode_length_slider_Callback(hObject, eventdata, handles)
length=get(hObject,'value');
set(handles.length_edit,'string',num2str(length));
guidata(hObject,handles);

function Electrode_length_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function prop_height_sld_Callback(hObject, eventdata, handles)
length=get(hObject,'value');
set(handles.prop_height_edit,'string',num2str(length));
guidata(hObject,handles);

function prop_height_sld_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function prop_height_edit_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.prop_height_sld,'value',str2num(edit));
len = str2num(get(handles.prop_height_edit,'String'));
minSliderValue = get(handles.prop_height_sld, 'Min');
maxSliderValue = get(handles.prop_height_sld, 'Max');
if len < minSliderValue
    l = minSliderValue;
    set(handles.prop_height_edit,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.prop_height_sld, 'Value', l);
elseif  len > maxSliderValue
     l = maxSliderValue;
    set(handles.prop_height_edit,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
   set(handles.prop_height_sld, 'Value', l);
end
guidata(hObject,handles);


function prop_height_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cap_inductance_sld_Callback(hObject, eventdata, handles)
i=get(hObject,'value');
set(handles.cap_inductance_edit,'string',num2str(i));
guidata(hObject,handles);

function cap_inductance_sld_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function cap_inductance_edit_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.cap_inductance_sld,'value',str2num(edit));
len = str2num(get(handles.cap_inductance_edit,'String'));
minSliderValue = get(handles.cap_inductance_sld, 'Min');
maxSliderValue = get(handles.cap_inductance_sld, 'Max');
if len < minSliderValue
    l = minSliderValue;
    set(handles.cap_inductance_edit,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.cap_inductance_sld, 'Value', l);
elseif  len > maxSliderValue
    l = maxSliderValue;
    set(handles.cap_inductance_edit,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.cap_inductance_sld, 'Value', l);
end
guidata(hObject,handles);

function cap_inductance_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function save_Callback(hObject, eventdata, handles)
state.Discharge_Energy_J=get(handles.Discharge_energy_edittext,'string');
state.Pulse_frequency_H=get(handles.Pulse_frequency_edit,'string');
state.Capacitance_uF=get(handles.Capacitance_edittext,'string');
state.Electrode_Inductance_nH=get(handles.electrode_inductance_edittext,'string');
state.Electrode_resistance_mOhm=get(handles.electrode_Resistance_edit,'string');
state.Capacitor_Inductance_nH=get(handles.cap_inductance_edit,'string');
state.Capacitor_resistance_mOhm=get(handles.ESR_edit,'string');

state.Electrode_spacing_mm=get(handles.ed_spacing_edit,'string');
state.Electrode_width_mm=get(handles.ed_width_edit,'string');
state.Electrode_length_mm=get(handles.length_edit,'string');
state.Electrode_thickness_mm=get(handles.electrode_thickness_edit,'string');
switch get(handles.Electrode_shape_popup,'Value')
    case 2
        state.Electrode_shape='Rectangular';
    case 3
       state.Electrode_shape='Tongue';
end

state.Activate_Flared_Geometry_micro_button=get(handles.Activate_Flared_Geometry_radiob,'value');
state.Flare_angle_degree=get(handles.Flare_angle_edittext,'string');
state.flared_length_mm=get(handles.flare_x_edit,'string');

state.Propellant_height_mm=get(handles.prop_height_edit,'string');
state.Propellant_width_mm=get(handles.Propellant_width_edittext,'string');
switch get(handles.Propellant_feeding_popup,'Value')
    case 2
       state.Propellant_feeding='Breech feeding';
    case 3
        state.Propellant_feeding='Side feeding';
    case 4
       state.Propellant_feeding='Combination feeding';
end

state.Spot_radius_mm=get(handles.spot_radius_micro,'string');
state.Spot_Velocity_km_s=get(handles.spot_velocity,'string');
state.Ion_erosion_rate_ug_C=get(handles.ion_erosion,'string');
state.Mean_ion_charge_state=get(handles.Mean_ion_charge,'string');
state.Ion_normalised_by_arc_current=get(handles.alpha_i,'string');
state.Current_per_spot_A=get(handles.current_per_spot,'string');

state.Exposed_propellant_area_cm2=get(handles.prop_area_edit,'string');
state.aspect_ratio=get(handles.ar_edit,'string');
state.EA_ratio_J_cm2=get(handles.EA_ratio_edittext,'string');
state.Electrode_erosion_per_Discharge_ug_s=get(handles.dm_dt,'string');

state.time_span_us=get(handles.time_results_edit,'string');
state.thrust_uN=get(handles.thrust_edit,'string');
state.Exhaust_velocity_m_s=get(handles.Exhaust_velocity_edit,'string');
state.Impulse_bit_uNs=get(handles.Impulse_edit,'string');
state.Specific_Impulse_s=get(handles.Specific_Impulse_edit,'string');
state.Inductance_gradient_nH_mm=get(handles.Inductance_gradient_edit,'string');
state.transfer_efficiency=get(handles.transfer_edit,'string');
state.acceleration_efficiency=get(handles.acceleration_edit,'string');
state.Thruster_efficiency=get(handles.Efficiency_edit,'string');
state.Inductance_change_nH=get(handles.Inductance_change_edit,'string');

startingFolder = 'C:\';
if ~exist(startingFolder, 'dir')
    startingFolder = pwd;% If that folder doesn't exist, just start in the current folder.
end
defaultFileName = fullfile(startingFolder, '*.mat');
[baseFileName, folder] = uiputfile(defaultFileName,'Select a mat file');% Put in the name of the mat file that the user wants to save, defaultFileName = fullfile(startingFolder, '*.mat')
if baseFileName == 0
    return;% User clicked the Cancel button.
end
fullFileName = fullfile(folder, baseFileName);
save(fullFileName,'state');


function insert_Callback(hObject, eventdata, handles)
startingFolder = 'C:\';
if ~exist(startingFolder, 'dir')
  startingFolder = pwd; % If that folder doesn't exist, just start in the current folder.
end
defaultFileName = fullfile(startingFolder, '*.mat');% Get the name of the mat file that the user wants to use.
[baseFileName, folder] = uigetfile(defaultFileName, 'Select a mat file');
if baseFileName == 0
  return; % User clicked the Cancel button. 
end
fullFileName = fullfile(folder, baseFileName);
load(fullFileName);
set(handles.Discharge_energy_edittext,'string',state.Discharge_Energy_J);
set(handles.Discharge_Energy_slider,'Value',str2num(state.Discharge_Energy_J));
set(handles.Pulse_frequency_edit,'string',state.Pulse_frequency_H);
set(handles.Pulse_frequency_slider,'Value',str2num(state.Pulse_frequency_H));
set(handles.Capacitance_edittext,'string',state.Capacitance_uF);
set(handles.Capacitance_slider,'Value',str2num(state.Capacitance_uF));
set(handles.electrode_inductance_edittext,'string',state.Electrode_Inductance_nH);
set(handles.electrode_Inductance_slider,'Value',str2num(state.Electrode_Inductance_nH));
set(handles.electrode_Resistance_edit,'string',state.Electrode_resistance_mOhm);
set(handles.electrode_Resistance_slider,'Value',str2num(state.Electrode_resistance_mOhm));
set(handles.cap_inductance_edit,'string',state.Capacitor_Inductance_nH);
set(handles.cap_inductance_sld,'Value',str2num(state.Capacitor_Inductance_nH));
set(handles.ESR_edit,'string',state.Capacitor_resistance_mOhm);
set(handles.ESR_slider,'Value',str2num(state.Capacitor_resistance_mOhm));

set(handles.ed_spacing_edit,'string',state.Electrode_spacing_mm);
set(handles.ed_spacing_slider,'Value',str2num(state.Electrode_spacing_mm));
set(handles.ed_width_edit,'string',state.Electrode_width_mm);
set(handles.ed_width_sld,'Value',str2num(state.Electrode_width_mm));
set(handles.length_edit,'string',state.Electrode_length_mm);
set(handles.Electrode_length_slider,'Value',str2num(state.Electrode_length_mm));
set(handles.electrode_thickness_edit,'string',state.Electrode_thickness_mm);
set(handles.electrode_thickness_sld,'Value',str2num(state.Electrode_thickness_mm));

if state.Electrode_shape=="Rectangular"
    set(handles.Electrode_shape_popup,'Value',2)
    else
    if state.Electrode_shape=="Tongue"
       set(handles.Electrode_shape_popup,'Value',3)
       
    end
end

button=state.Activate_Flared_Geometry_micro_button;
if button==0
    set(handles.Activate_Flared_Geometry_radiob,'value',0)
    handlesArr = [handles.Flare_Angle_slider,handles.Flare_angle_edittext,handles.flare_x_edit,handles.flare_x_sld, handles.flare_angle_slider_large, handles.flare_angle_edit_large, handles.flare_length_slider_large, handles.flare_length_edit_large];
    set(handlesArr, 'Enable', 'off');
    set(handles.Flare_angle_edittext,'string','null');
    set(handles.Flare_Angle_slider,'Value',str2num('0'));
    set(handles.flare_x_edit,'string','null');
    set(handles.flare_x_sld,'Value',str2num('0'));
else
    set(handles.Activate_Flared_Geometry_radiob,'value',1)
    handlesArr = [handles.Flare_Angle_slider,handles.Flare_angle_edittext,handles.flare_x_edit,handles.flare_x_sld, handles.flare_angle_slider_large, handles.flare_angle_edit_large, handles.flare_length_slider_large, handles.flare_length_edit_large];
    set(handlesArr, 'Enable', 'on');
    set(handles.Flare_angle_edittext,'string',state.Flare_angle_degree);
    set(handles.Flare_Angle_slider,'Value',str2num(state.Flare_angle_degree));
    set(handles.flare_x_edit,'string',state.flared_length_mm);
    set(handles.flare_x_sld,'Value',str2num(state.flared_length_mm));
end
set(handles.prop_height_edit,'string',state.Propellant_height_mm);
set(handles.prop_height_sld,'Value',str2num(state.Propellant_height_mm));
set(handles.Propellant_width_edittext,'string',state.Propellant_width_mm);
set(handles.Propellant_width_slider,'Value',str2num(state.Propellant_width_mm));

if state.Propellant_feeding=="Breech feeding"
    set(handles.Propellant_feeding_popup,'Value',2)
else 
    if state.Propellant_feeding=="Side feeding"
        set(handles.Propellant_feeding_popup,'Value',3)
    else
        if state.Propellant_feeding=="Combination feeding"
            set(handles.Propellant_feeding_popup,'Value',4)
        end     
    end
end

set(handles.spot_radius_micro,'string',state.Spot_radius_mm);
set(handles.spot_velocity,'string',state.Spot_Velocity_km_s);
set(handles.ion_erosion,'string',state.Ion_erosion_rate_ug_C);
set(handles.Mean_ion_charge,'string',state.Mean_ion_charge_state);
set(handles.alpha_i,'string',state.Ion_normalised_by_arc_current);
set(handles.current_per_spot,'string',state.Current_per_spot_A);

set(handles.prop_area_edit,'string',state.Exposed_propellant_area_cm2);
set(handles.ar_edit,'string',state.aspect_ratio);
set(handles.EA_ratio_edittext,'string',state.EA_ratio_J_cm2);
set(handles.dm_dt,'string',state.Electrode_erosion_per_Discharge_ug_s);

set(handles.time_results_edit,'string',state.time_span_us);
set(handles.thrust_edit,'string',state.thrust_uN);
set(handles.Exhaust_velocity_edit,'string',state.Exhaust_velocity_m_s);
set(handles.Impulse_edit,'string',state.Impulse_bit_uNs);
set(handles.Specific_Impulse_edit,'string',state.Specific_Impulse_s);
set(handles.Inductance_gradient_edit,'string',state.Inductance_gradient_nH_mm);
set(handles.transfer_edit,'string',state.transfer_efficiency);
set(handles.acceleration_edit,'string',state.acceleration_efficiency);
set(handles.Efficiency_edit,'string',state.Thruster_efficiency);
set(handles.Inductance_change_edit,'string',state.Inductance_change_nH);


function plot_choice_Callback(hObject, eventdata, handles)

function plot_choice_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function figure1_SizeChangedFcn(hObject, eventdata, handles)

function main_menu_Callback(hObject, eventdata, handles)
set(handles.main_panel,'Visible','On');
set(handles.main_menu,'BackgroundColor',[1,1,1]);
set(handles.micro_ppt,'BackgroundColor',[.8,.8,.8]);
set(handles.large_ppt,'BackgroundColor',[.8,.8,.8]);
set(handles.micro_ppt_panel,'Visible','Off');
set(handles.large_ppt_panel,'Visible','Off');

function micro_ppt_Callback(hObject, eventdata, handles)
set(handles.micro_ppt,'BackgroundColor',[1,1,1]);
set(handles.large_ppt,'BackgroundColor',[.8,.8,.8]);
set(handles.main_menu,'BackgroundColor',[.8,.8,.8]);
set(handles.main_panel,'Visible','Off');
set(handles.micro_ppt_panel,'Visible','On');
set(handles.large_ppt_panel,'Visible','Off');

function large_ppt_Callback(hObject, eventdata, handles)
set(handles.large_ppt_panel,'Visible','On');
set(handles.large_ppt,'BackgroundColor',[1,1,1]);
set(handles.micro_ppt,'BackgroundColor',[.8,.8,.8]);
set(handles.main_menu,'BackgroundColor',[.8,.8,.8]);
set(handles.main_panel,'Visible','Off');
set(handles.micro_ppt_panel,'Visible','Off');


function pulse_frequency_slider_large_Callback(hObject, eventdata, handles)
Pfreq=get(hObject,'value');
set(handles.pulse_frequency_edit_large,'string',num2str(Pfreq));
guidata(hObject,handles);

function pulse_frequency_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function capacitance_slider_large_Callback(hObject, eventdata, handles)
Cap=get(hObject,'value');
set(handles.capacitance_edit_large,'string',num2str(Cap));
guidata(hObject,handles);

function capacitance_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function discharge_energy_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.discharge_energy_slider_large,'value',str2num(edit));
energy = str2num(get(handles.discharge_energy_edit_large,'String'));
minSliderValue = get(handles.discharge_energy_slider_large, 'Min');
maxSliderValue = get(handles.discharge_energy_slider_large, 'Max');
if energy < minSliderValue
     e = minSliderValue ;
    set(handles.discharge_energy_edit_large,'String',num2str(e));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.discharge_energy_slider_large, 'Value', e);
elseif  energy > maxSliderValue
    e = maxSliderValue;
    set(handles.discharge_energy_edit_large,'String',num2str(e));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.discharge_energy_slider_large, 'Value', e);
end
guidata(hObject,handles);

function discharge_energy_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function discharge_energy_slider_large_Callback(hObject, eventdata, handles)
Den=get(hObject,'value');
set(handles.discharge_energy_edit_large,'string',num2str(Den));
guidata(hObject,handles);

function discharge_energy_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function pulse_frequency_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.pulse_frequency_slider_large,'value',str2num(edit));
freq = str2num(get(handles.pulse_frequency_edit_large,'String'));
minSliderValue = get(handles.pulse_frequency_slider_large, 'Min');
maxSliderValue = get(handles.pulse_frequency_slider_large, 'Max');
if freq < minSliderValue
    f = minSliderValue;
    set(handles.pulse_frequency_edit_large,'String',num2str(f));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.pulse_frequency_slider_large, 'Value', f);
elseif  freq > maxSliderValue 
    f = maxSliderValue;
    set(handles.pulse_frequency_edit_large,'String',num2str(f));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.pulse_frequency_slider_large, 'Value', f);
end
guidata(hObject,handles);

function pulse_frequency_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function save_button_large_Callback(hObject, eventdata, handles)
state.Discharge_energy_J=get(handles.discharge_energy_edit_large,'string');
state.Pulse_frequency_Hz=get(handles.pulse_frequency_edit_large,'string');
state.Capacitance_uF=get(handles.capacitance_edit_large,'string');
state.Electrode_Inductance_nH=get(handles.electrode_inductance_edit_large,'string');
state.Electrode_resistance_mOhm=get(handles.electrode_resistance_edit_large,'string');
state.Capcitor_Inductance_nH=get(handles.capacitor_inductance_edit_large,'string');
state.Capacitor_resistance_mOhm=get(handles.capacitor_resistance_edit_large,'string');

state.Electrode_spacing_mm=get(handles.electrode_spacing_edit_large,'string');
state.Electrde_width_mm=get(handles.electrode_width_edit_large,'string');
state.Electrode_length_mm=get(handles.electrode_length_edit_large,'string');
state.electrode_thickness_mm=get(handles.electrode_thickness_edit_large,'string');

state.Propellant_height_mm=get(handles.propellant_height_edit_large,'string');
state.Propellant_width_mm=get(handles.propellant_width_edit_large,'string');
switch get(handles.electrode_shape_popup_large,'Value')
    case 2
        state.Electrode_shape='Rectangular';
    case 3
       state.Electrode_shape='Tongue';
end

switch get(handles.propellant_feed_popup_large,'Value')
    case 2
       state.Propellant_feeding='Breech feeding';
    case 3
        state.Propellant_feeding='Side feeding';
    case 4
       state.Propellant_feeding='Combination feeding';
end

state.Activate_Flared_Geometry_large_button=get(handles.flare_button_large,'value');
state.Flared_angle_degree=get(handles.flare_angle_edit_large,'string');
state.Flared_length_mm=get(handles.flare_length_edit_large,'string');

state.Spot_radius_mm=get(handles.spot_radius_large,'string');
state.Spot_Velocity_km_s=get(handles.spot_velocity_large,'string');
state.Ion_erosion_rate_ug_C=get(handles.ion_erosion_large,'string');
state.Mean_ion_charge_state=get(handles.mean_ion_charge_large,'string');
state.Ion_normalised_by_arc_current=get(handles.alpha_large,'string');
state.Current_per_spot_A=get(handles.spot_current_large,'string');

state.Time_results_us=get(handles.time_edit_large,'string');
state.Exposed_Propellant_area_cm2=get(handles.exposed_propellant_large,'string');
state.Aspect_ratio=get(handles.aspect_ratio_large,'string');
state.E_A_ratio_J_cm2=get(handles.E_A_ratio_large,'string');
state.Electrode_erosion_per_Discharge_ug_s=get(handles.electrode_erosion_large,'string');
state.Thrust_uN=get(handles.thrust_large,'string');
state.Exhaust_velocity_m_s=get(handles.exhaust_velocity_large,'string');
state.Impulse_bit_uNs=get(handles.impulse_bit_large,'string');
state.Specific_Impulse_s=get(handles.specific_impulse_large,'string');
state.Inductance_gradient_uH_mm=get(handles.inductance_gradient_large,'string');
state.Transfer_Efficiency=get(handles.transfer_efficiency_large,'string');
state.Acceleration_Efficiency=get(handles.acceleration_eff_large,'string');
state.Efficiency=get(handles.thruster_eff_large,'string');
state.Inductance_change_nH=get(handles.inductance_change_large,'string');

startingFolder = 'C:\';
if ~exist(startingFolder, 'dir')
    startingFolder = pwd;% If that folder doesn't exist, just start in the current folder.
end
defaultFileName = fullfile(startingFolder, '*.mat');
[baseFileName, folder] = uiputfile(defaultFileName,'Select a mat file');% Put in the name of the mat file that the user wants to save, defaultFileName = fullfile(startingFolder, '*.mat')
if baseFileName == 0
    return; % User clicked the Cancel button.
end
fullFileName = fullfile(folder, baseFileName);
save(fullFileName,'state');

function insert_large_Callback(hObject, eventdata, handles)
startingFolder = 'C:\';
if ~exist(startingFolder, 'dir')
  startingFolder = pwd;% If that folder doesn't exist, just start in the current folder.
end
defaultFileName = fullfile(startingFolder, '*.mat');% Get the name of the mat file that the user wants to use.
[baseFileName, folder] = uigetfile(defaultFileName, 'Select a mat file');
if baseFileName == 0
  return;% User clicked the Cancel button. 
end
fullFileName = fullfile(folder, baseFileName);
load(fullFileName);
set(handles.discharge_energy_edit_large,'string',state.Discharge_energy_J);
set(handles.discharge_energy_slider_large,'value',str2num(state.Discharge_energy_J));
set(handles.pulse_frequency_edit_large,'string',state.Pulse_frequency_Hz);
set(handles.pulse_frequency_slider_large,'value',str2num(state.Pulse_frequency_Hz));
set(handles.capacitance_edit_large,'string',state.Capacitance_uF);
set(handles.capacitance_slider_large,'value',str2num(state.Capacitance_uF));
set(handles.electrode_inductance_edit_large,'string',state.Electrode_Inductance_nH);
set(handles.electrode_inductance_slider_large,'value',str2num(state.Electrode_Inductance_nH));
set(handles.electrode_resistance_edit_large,'string',state.Electrode_resistance_mOhm);
set(handles.electrode_resistance_slider_large,'value',str2num(state.Electrode_resistance_mOhm));
set(handles.capacitor_inductance_edit_large,'string',state.Capcitor_Inductance_nH);
set(handles.capacitor_inductance_slider_large,'value',str2num(state.Capcitor_Inductance_nH));
set(handles.capacitor_resistance_edit_large,'string',state.Capacitor_resistance_mOhm);
set(handles.capacitor_resiatnce_slider_large,'value',str2num(state.Capacitor_resistance_mOhm));

set(handles.electrode_spacing_edit_large,'string',state.Electrode_spacing_mm);
set(handles.electrode_spacing_slider_large,'value',str2num(state.Electrode_spacing_mm));
set(handles.electrode_width_edit_large,'string',state.Electrde_width_mm);
set(handles.electrode_width_slider_large,'value',str2num(state.Electrde_width_mm));
set(handles.electrode_length_edit_large,'string',state.Electrode_length_mm);
set(handles.electrode_length_slider_large,'value',str2num(state.Electrode_length_mm));
set(handles.electrode_thickness_edit_large,'string',state.electrode_thickness_mm);
set(handles.electrode_thickness_slider_large,'value',str2num(state.electrode_thickness_mm));

set(handles.propellant_height_edit_large,'string',state.Propellant_height_mm);
set(handles.propellant_height_slider_large,'value',str2num(state.Propellant_height_mm));
set(handles.propellant_width_edit_large,'string',state.Propellant_width_mm);
set(handles.propellant_width_slider_large,'value',str2num(state.Propellant_width_mm));

if state.Electrode_shape=="Rectangular"
     ;set(handles.electrode_shape_popup_large,'Value',2);
else
    if state.Electrode_shape=="Tongue"
        set(handles.electrode_shape_popup_large,'Value',3);
    end
end

if state.Propellant_feeding=="Breech feeding"
    set(handles.propellant_feed_popup_large,'Value',2);  
else
    if state.Propellant_feeding=="Side feeding"
        set(handles.propellant_feed_popup_large,'Value',3);
    else
        if state.Propellant_feeding=="Combination feeding"
            set(handles.propellant_feed_popup_large,'Value',4);
        end
    end
end

button=state.Activate_Flared_Geometry_large_button;
if button==0
    set(handles.flare_button_large,'value',0);
    handlesArr=[handles.flare_angle_slider_large, handles.flare_angle_edit_large,handles.flare_length_slider_large, handles.flare_length_edit_large];
    set(handlesArr, 'Enable', 'off');
    set(handles.flare_angle_edit_large,'string','null');
    set(handles.flare_angle_slider_large,'value',str2num('0'));
    set(handles.flare_length_edit_large,'string','null');
    set(handles.flare_length_slider_large,'value',str2num('0'));
else
    set(handles.flare_button_large,'value',1);
    handlesArr=[handles.flare_angle_slider_large, handles.flare_angle_edit_large,handles.flare_length_slider_large, handles.flare_length_edit_large];
    set(handlesArr, 'Enable', 'on');
    set(handles.flare_angle_edit_large,'string',state.Flared_angle_degree);
    set(handles.flare_angle_slider_large,'value',str2num(state.Flared_angle_degree));
    set(handles.flare_length_edit_large,'string',state.Flared_length_mm);
    set(handles.flare_length_slider_large,'value',str2num(state.Flared_length_mm));
end
set(handles.spot_radius_large,'string',state.Spot_radius_mm);
set(handles.spot_velocity_large,'string',state.Spot_Velocity_km_s);
set(handles.ion_erosion_large,'string',state.Ion_erosion_rate_ug_C);
set(handles.mean_ion_charge_large,'string',state.Mean_ion_charge_state);
set(handles.alpha_large,'string',state.Ion_normalised_by_arc_current);
set(handles.spot_current_large,'string',state.Current_per_spot_A);

set(handles.time_edit_large,'string',state.Time_results_us);
set(handles.exposed_propellant_large,'string',state.Exposed_Propellant_area_cm2);
set(handles.aspect_ratio_large,'string',state.Aspect_ratio);
set(handles.E_A_ratio_large,'string',state.E_A_ratio_J_cm2);
set(handles.electrode_erosion_large,'string',state.Electrode_erosion_per_Discharge_ug_s);
set(handles.thrust_large,'string',state.Thrust_uN);
set(handles.exhaust_velocity_large,'string',state.Exhaust_velocity_m_s);
set(handles.impulse_bit_large,'string',state.Impulse_bit_uNs);
set(handles.specific_impulse_large,'string',state.Specific_Impulse_s);
set(handles.inductance_gradient_large,'string',state.Inductance_gradient_uH_mm);
set(handles.transfer_efficiency_large,'string',state.Transfer_Efficiency);
set(handles.acceleration_eff_large,'string',state.Acceleration_Efficiency);
set(handles.thruster_eff_large,'string',state.Efficiency);
set(handles.inductance_change_large,'string',state.Inductance_change_nH);


function reset_button_large_Callback(hObject, eventdata, handles)
handles1=[handles.capacitor_resiatnce_slider_large,handles.electrode_thickness_slider_large,handles.flare_angle_slider_large,handles.propellant_height_slider_large,handles.capacitor_resistance_edit_large,handles.discharge_energy_slider_large,handles.pulse_frequency_slider_large,handles.capacitance_slider_large,handles.electrode_inductance_slider_large,handles.electrode_resistance_slider_large,handles.capacitor_inductance_slider_large];
set(handles1,'value',0);
handles2=[handles.range_min_large,handles.range_max_large,handles.capacitor_resistance_edit_large,handles.propellant_height_edit_large,handles.discharge_energy_edit_large,handles.pulse_frequency_edit_large,handles.capacitance_edit_large,handles.electrode_inductance_edit_large,handles.electrode_resistance_edit_large,handles.capacitor_inductance_edit_large];
set(handles2,'string','');
handles3=[handles.flare_angle_edit_large,handles.electrode_spacing_edit_large,handles.electrode_width_edit_large,handles.electrode_length_edit_large,handles.flare_length_edit_large,handles.exposed_propellant_large,handles.aspect_ratio_large,handles.E_A_ratio_large,handles.electrode_thickness_edit_large];
set(handles3,'string','');
handles4=[handles.resolution_large,handles.time_edit_large,handles.thrust_large,handles.exhaust_velocity_large,handles.impulse_bit_large,handles.specific_impulse_large,handles.inductance_gradient_large,handles.transfer_efficiency_large,handles.acceleration_eff_large,handles.thruster_eff_large,handles.inductance_change_large];
set(handles4,'string','');
handles5=[handles.propellant_width_slider_large,handles.electrode_spacing_slider_large,handles.electrode_width_slider_large,handles.electrode_length_slider_large,handles.flare_length_slider_large];
set(handles5,'value',0);
handles6=[handles.electrode_shape_popup_large,handles.dependent_variable_popup_large,handles.independent_variable_popup_large,handles.propellant_feed_popup_large,handles.plot_popup_large];
set(handles6,'value',1);
handles7=[handles.spot_radius_large,handles.spot_velocity_large,handles.ion_erosion_large,handles.mean_ion_charge_large,handles.alpha_large,handles.spot_current_large,handles.propellant_width_edit_large,handles.electrode_erosion_large];
set(handles7,'string','');
axes(handles.axes1_large);
cla reset;
axes(handles.axes2_large);
cla reset;
axes(handles.axes3_large);
cla reset;
set(handles.flare_button_large,'value',0);

function thrust_large_Callback(hObject, eventdata, handles)

function thrust_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function exhaust_velocity_large_Callback(hObject, eventdata, handles)

function exhaust_velocity_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function impulse_bit_large_Callback(hObject, eventdata, handles)

function impulse_bit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function thruster_eff_large_Callback(hObject, eventdata, handles)

function thruster_eff_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function inductance_gradient_large_Callback(hObject, eventdata, handles)

function inductance_gradient_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function specific_impulse_large_Callback(hObject, eventdata, handles)

function specific_impulse_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function transfer_efficiency_large_Callback(hObject, eventdata, handles)

function transfer_efficiency_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function acceleration_eff_large_Callback(hObject, eventdata, handles)

function acceleration_eff_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function inductance_change_large_Callback(hObject, eventdata, handles)

function inductance_change_large_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function clear_ax1_large_Callback(hObject, eventdata, handles)
set(handles.clear_ax1_large,'BackgroundColor',[1,1,1]);
axes(handles.axes1_large);
cla reset;
set(handles.clear_ax1_large,'BackgroundColor',[0.85,0.85,0.85]);

function save_ax1_large_Callback(hObject, eventdata, handles)
set(handles.save_ax1_large,'BackgroundColor',[1,1,1]);
set(handles.axes1_large,'units','normalized');
fignew = figure('Visible','on'); % Invisible figure
newAxes = copyobj(handles.axes1_large,fignew); % Copy the appropriate axes
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
set(fignew,'CreateFcn','set(gcbf,''Visible'',''on'')'); % Make it visible upon loading
set(fignew, 'PaperPositionMode', 'auto');  
set(handles.save_ax1_large,'BackgroundColor',[0.85,0.85,0.85]);


function clear_ax2_large_Callback(hObject, eventdata, handles)
set(handles.clear_ax2_large,'BackgroundColor',[1,1,1]);
axes(handles.axes2_large);
cla reset;
set(handles.clear_ax2_large,'BackgroundColor',[0.85,0.85,0.85]);

function save_ax2_large_Callback(hObject, eventdata, handles)
set(handles.save_ax2_large,'BackgroundColor',[1,1,1]);
set(handles.axes2_large,'units','normalized');
fignew = figure('Visible','on'); % Invisible figure
newAxes = copyobj(handles.axes2_large,fignew); % Copy the appropriate axes
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
set(fignew,'CreateFcn','set(gcbf,''Visible'',''on'')'); % Make it visible upon loading
set(fignew, 'PaperPositionMode', 'auto'); 
set(handles.save_ax2_large,'BackgroundColor',[0.85,0.85,0.85]);

function clear_ax3_large_Callback(hObject, eventdata, handles)
set(handles.clear_ax3_large,'BackgroundColor',[1,1,1]);
axes(handles.axes3_large);
cla reset;
set(handles.clear_ax3_large,'BackgroundColor',[0.85,0.85,0.85]);

function save_ax3_large_Callback(hObject, eventdata, handles)
set(handles.save_ax3_large,'BackgroundColor',[1,1,1]);
set(handles.axes3_large,'units','normalized');
fignew = figure('Visible','on'); % Invisible figure
newAxes = copyobj(handles.axes3_large,fignew); % Copy the appropriate axes
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
set(fignew,'CreateFcn','set(gcbf,''Visible'',''on'')'); % Make it visible upon loading
set(fignew, 'PaperPositionMode', 'auto'); 
set(handles.save_ax3_large,'BackgroundColor',[0.85,0.85,0.85]);


function plot_button_large_Callback(hObject, eventdata, handles)
set(handles.plot_button_large,'BackgroundColor',[1,1,1]);
digits(7);
set(handles.large_ppt,'BackgroundColor',[1,1,1]);
progress_f= waitbar(0,'Please wait...');
waitbar(.33,progress_f,'Loading your data');
pause(0.5);
switch get(handles.plot_popup_large,'Value')
    case 2
        axes(handles.axes1_large);
    case 3
        axes(handles.axes2_large);
    case 4
        axes(handles.axes3_large);
end
min=(str2num(get(handles.range_min_large,'string')))*1e-3;
max=(str2num(get(handles.range_max_large,'string')))*1e-3;
res=(str2num(get(handles.resolution_large,'string')))*1e-3;
%calculations
u=1.2566370614*(1e-6);%permeability of free space
g=9.8;
E=str2num(get(handles.discharge_energy_edit_large,'string'));%discharge energy in J
f=str2num(get(handles.pulse_frequency_edit_large,'string'));%pulse frequency in Hz
C=(str2num(get(handles.capacitance_edit_large,'string')))*1e-6;%capacitance in F
Le=(str2num(get(handles.electrode_inductance_edit_large,'string')))*1e-9;%electrode inductance in H
Re=(str2num(get(handles.electrode_resistance_edit_large,'string')))*1e-3;%electrode resistance in ohm                              
ESR=(str2num(get(handles.capacitor_resistance_edit_large,'string')))*1e-3;%capacitor resistance in ohm 
Lcap=(str2num(get(handles.capacitor_inductance_edit_large,'string')))*1e-9;%capacitor inductance in H
h=(str2num(get(handles.electrode_spacing_edit_large,'string')))*1e-3;%electrode spacing
w=(str2num(get(handles.electrode_width_edit_large,'string')))*1e-3;%electrode width
l=(str2num(get(handles.electrode_length_edit_large,'string')))*1e-3;%electrode length
d=(str2num(get(handles.electrode_thickness_edit_large,'string')))*1e-3;%electrode thickness
t_span=(str2num(get(handles.time_edit_large,'string')))*1e-6;%time_span for plasma sheet to reach the end of electrodes                                
if (h/w <2)
    warndlg('It seems the aspect ratio is less than 1. This makes your calculation into the catrgory of Micro PPT. Kindly, check the dimensions again or use the Micro PPT caculation page.');
else  
    L_large= 0.4*l*(log(h/(w+d)) + 3/2 - h/l +0.22*(w+d)/l);
    L=(L_large*1e-6)+Lcap+Le;%total initial resistance
    R=ESR+Re;%total initial resistance
    Vo=sqrt(2*E/C);%voltage
    zi=(R/2)*sqrt(C/L);%zeta
    wo=1/sqrt(L*C);
    syms t;%declaring t as a variaable
    if zi<1
        i2=@(t)((Vo^2.*exp(-2*zi*wo*t)*(sin(sqrt(1-zi^2)*wo*t))^2)/(L^2.*(1-zi^2).*wo^2));%current equation
        i=@(t)((Vo.*exp(-zi*wo*t)*(sin(sqrt(1-(zi^2)).*wo*t)))/(L.*sqrt(1-(zi^2)).*wo));
        v=@(t)((Vo.*zi.*exp(-zi*wo*t).*sin(sqrt(1-(zi^2)).*wo*t)./sqrt(1-(zi^2)))+(Vo.*exp(-zi*wo*t).*cos(sqrt(1-(zi^2)).*wo*t))); %voltage equation
    else
        if C==(4*L/(R^2))
            i2=@(t)(Vo^2/L^2)*(t^2.*exp(-2*wo*t));%current equation
            i=@(t)(Vo/L)*(t.*exp(-wo*t));
            v=@(t)(Vo.*exp(-wo*t).*(1+wo*t));%voltage equation
        else
            i2=@(t)(Vo^2.*(exp(-2*zi*wo*t))*(sinh(wo*t*sqrt(zi^2-1)))^2)/(L^2.*(zi^2-1).*wo^2);%current equation
            i=@(t)(Vo.*(exp(-zi*wo*t))*(sinh(wo*t*sqrt((zi^2)-1))))/(L.*sqrt((zi^2)-1).*wo);
            v=@(t)((Vo.*zi.*exp(-zi*wo*t).*sinh(sqrt((zi^2)-1).*wo*t)./sqrt((zi^2)-1))+(Vo.*exp(-zi*wo*t).*cosh(sqrt((zi^2)-1).*wo*t)));
        end
    end
    V=v(t_span);
    I=i(t_span);%loop current flowing through the circuit
    I_arc=Vo*sqrt(C/L);%arc current flowing through the circuit loop
    i2_integration=int(i2,t,[0,t_span]);%evaluation of integral of (i(t))^2
    Xp=vpa(i2_integration);%making the above answer definite
    %calculation of mbit
    me=9.1.*1e-31;%mass of electron
    e=1.6.*1e-19;%elementary charge
    eo=8.85.*1e-12;%Vacuum permittivity of free space(F/m)
    Ro=(1/3)*sqrt(C/L);
    Lp=0;%plasma inductance
    Te=11605;%temperature(in K)(1 eV)
    r=(str2num(get(handles.spot_radius_large,'string')))*1e-3;%radius of the cathode spot
    S_spot=pi.*r^2;%Surface area of the cathode spot near the mixing region
    I_spot=(str2num(get(handles.spot_current_large,'string')));%For copper electrodes it has been experimentally observed that the current per observed spot
    So=S_spot*I_arc/I_spot;%Surface area of the initial plasma flow near the mixing region
    del_T= 8*11605;%temperature difference(8 eV)(in K)
    lambda=1.2*del_T;%latent heat of vaporization, specific heat=1.2
    alp_I=(str2num(get(handles.alpha_large,'string')));%Ion current normalised by arc current, the ratio between the ion current and the arc current in a cathodic plasma
    Q=(str2num(get(handles.mean_ion_charge_large,'string')));%mean ion charge state number in plasma flow
    V_spot=(str2num(get(handles.spot_velocity_large,'string')))*1e3;%initial copper plasma has a velocity towards the direction of the anode
    Ne=I*alp_I./(Q*e*So*V_spot);%Electron density
    Ni=Ne/Q;%ion density
    A=23-log(Ne^0.5*Q*Te^-1.5);%Coulomb logarithm, factor by which small angle collisions are more effective than large angle collision
    neu=(3.62.*1e-6)*A*Ne*(Te^(-3/2))/Q;%The frequency at which electrons and ions collide
    sigma=(Ne.*(e^2))/(me*neu);%conductivity
    Rp=abs(V/I);%plasma resistance
    V_sheath=abs(V/2);
    d_sheath=sqrt(2*V_sheath*eo/(e*Ni));          %thickness of the sheet
    waitbar(.67,progress_f,'Processing your data');
    pause(0.5);
    syms x;
    switch get(handles.independent_variable_popup_large,'Value')
        case 2
            h = min:res:max;
            %for non flared geometry
            if (get(handles.flare_button_large,'value'))==0
                m_bit=(sigma*d_sheath*d./(lambda.*h)).*(Rp^2+Lp^2/L*C)*(E/Ro);
                mbit=vpa(m_bit);
                switch get(handles.electrode_shape_popup_large,'Value')
                    case 2
                        xp=(u.*h./(2.*mbit.*w)).*Xp;
                        ind_grad=(0.6+0.4.*log(h./(w+d)))*1e-6;
                    case 3
                        tl=0.2*l;%part with electrode in tongue shape
                        xp1=(u.*h./(2.*mbit.*w)).*Xp;%calculation for rectangular part
                        xp1=vpa(xp1);
                        xp2=@(x)(((u.*h./(2.*mbit.*w)).*Xp)./(1-x/l));%calculation for tongue part
                        xp2=vpa(int(xp2,x,[0,tl]));
                        xp=xp1+xp2;
                        ind_grad=(0.6+0.4.*log(h./(w+d)))*1e-6;
                end
            end
            %for flared geometry
            if (get(handles.flare_button_large,'value'))==1
                flare_l=((str2num(get(handles.flare_length_edit_large,'string')))*1e-3);
                alpha= str2num(get(handles.flare_angle_edit_large,'string'));%flare angle
                switch get(handles.electrode_shape_popup_large,'Value')
                    case 2
                        m_bit1=(sigma*d_sheath*d./(lambda.*h)).*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for rectangular parallel part
                        m_bit1=vpa(m_bit1);
                        m_bit2=@(x)((sigma*d_sheath*d./(lambda.*(h+(tand(alpha))*x))).*(Rp^2+Lp^2./(L*C)).*(E/Ro));%calculation for rectangular flared part
                        m_bit2=vpa(int(m_bit2,x,[0,flare_l]));
                        mbit=m_bit1+m_bit2;
                        xp1=(u.*h./(2.*m_bit1.*w)).*Xp;%calculation for rectangular parallel part
                        xp1=vpa(xp1);
                        xp2=@(x)(((u.*(h+(tand(alpha))*x))./(2.*m_bit2*w)).*Xp);%calculation for rectangular flared part
                        xp2=vpa(int(xp2,x,[0,flare_l]));
                        xp=xp1+xp2;
                        ind_grad1=(0.6+0.4.*log(h./(w+d)))*1e-6;%calculation for rectangular parallel part
                        ind_grad2=@(x)(0.6+0.4.*log((h+tand(alpha)*x)./(w+d)))*1e-6;
                        ind_grad2=vpa(int(ind_grad2,x,[0,flare_l]));%calculation for rectangular flared part
                        ind_grad=ind_grad1+ind_grad2;
                    case 3
                        m_bit1=(sigma*d_sheath*d./(lambda.*h))*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for rectangular parallel part
                        m_bit1=vpa(m_bit1);
                        m_bit2=@(x)((sigma*d_sheath.*d)./(lambda.*(h+(tand(alpha))*x))).*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for flared part
                        m_bit2=vpa(int(m_bit2,x,[0,flare_l]));
                        mbit=m_bit1+m_bit2;
                        xp1=(u.*h./(2.*m_bit1.*w)).*Xp;%calculation for rectangular parallel part
                        flare_l_a=0.8.*flare_l;%part with electrode in rectangular shape
                        flare_l_b=0.2.*flare_l;%part with electrode in tongue shape
                        xp2=@(x)(((u.*(h+(tand(alpha))*x))./(2.*m_bit2.*w)).*Xp);
                        xp2=vpa(int(xp2,x,[0,flare_l_a]));
                        xp3=@(x)(((u.*(h+(tand(alpha))*x))./(2.*m_bit2.*(w*(1-x/l)))).*Xp);
                        xp3=vpa(int(xp3,x,[0,flare_l_b]));
                        xp=xp1+xp2+xp3;
                        ind_grad1=(0.6+0.4.*log(h./(w+d)))*1e-6;%calculation for rectangular parallel part
                        ind_grad2=@(x)(0.6+0.4.*log((h+tand(alpha)*x)./(w+d)))*1e-6;%calculation for rectangular flared part
                        ind_grad2=vpa(int(ind_grad2,x,[0,flare_l_a]));
                        ind_grad3=@(x)(0.6+0.4.*log((h+tand(alpha)*x)./((w.*(1-x/l))+d)))*1e-6;%calculation for tongue flared part
                        ind_grad3=vpa(int(ind_grad3,x,[0,flare_l_b]));
                        ind_grad=ind_grad1+ind_grad2+ind_grad3;
                end
            end
            ve=xp;%exhaust velocity
            Isp = ve/g;
            Y=1.3;%for teflon
            impulse_bit_EM=Isp.*mbit.*g.*f; %impulse bit
            impulse_bit_ET=sqrt((8*(Y-1).*mbit.*E)/(Y^2.*(Y+1)))+(Xp.*ind_grad)/2;%impulse bit
            impulse_bit= impulse_bit_EM+ impulse_bit_ET;
            Efficiency=impulse_bit.*ve./(2*E); %efficiency
            thrust=impulse_bit*f;
            h=h*1e3;
            switch get(handles.dependent_variable_popup_large,'Value')
                case 2
                    plot(h,Efficiency*1e2,'-*');
                    xlabel('electrode spacing(mm)');
                    ylabel('Thruster Efficiency(%)');
                case 3
                    plot(h,impulse_bit*1e6,'-*');
                    xlabel('electrode spacing(mm)');
                    ylabel('Impulse Bit(uNs)');
                case 4
                    plot(h,ve,'-*');
                    xlabel('electrode spacing(mm)');
                    ylabel('Exhaust Velocity(m/s)');
                case 5
                    plot(h,thrust*1e6,'-*');
                    xlabel('electrode spacing(mm)');
                    ylabel('Thrust(uN)');
            end
        case 3
            w = min:res:max;
            %for non flared geometry
            if (get(handles.flare_button_large,'value'))==0
                m_bit=(sigma*d_sheath*d./(lambda.*h)).*(Rp^2+Lp^2/L*C)*(E/Ro);
                mbit=vpa(m_bit);
                switch get(handles.electrode_shape_popup_large,'Value')
                    case 2
                        xp=(u.*h./(2.*mbit.*w)).*Xp;
                        ind_grad=(0.6+0.4.*log(h./(w+d)))*1e-6;
                    case 3
                        tl=0.2*l;%part with electrode in tongue shape
                        xp1=(u.*h./(2.*mbit.*w)).*Xp;%calculation for rectangular part
                        xp1=vpa(xp1);
                        xp2=@(x)(((u.*h./(2.*mbit.*w)).*Xp)./(1-x/l));%calculation for tongue part
                        xp2=vpa(int(xp2,x,[0,tl]));
                        xp=xp1+xp2;
                        ind_grad=(0.6+0.4.*log(h./(w+d)))*1e-6;
                end
            end
            %for flared geometry
            if (get(handles.flare_button_large,'value'))==1
                flare_l=((str2num(get(handles.flare_length_edit_large,'string')))*1e-3);
                alpha=str2num(get(handles.flare_angle_edit_large,'string'));%flare angle
                switch get(handles.electrode_shape_popup_large,'Value')
                    case 2
                        m_bit1=(sigma*d_sheath*d./(lambda.*h)).*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for rectangular parallel part
                        m_bit1=vpa(m_bit1);
                        m_bit2=@(x)((sigma*d_sheath*d./(lambda.*(h+(tand(alpha))*x))).*(Rp^2+Lp^2./(L*C)).*(E/Ro));%calculation for rectangular flared part
                        m_bit2=vpa(int(m_bit2,x,[0,flare_l]));
                        mbit=m_bit1+m_bit2;
                        xp1=(u.*h./(2.*m_bit1.*w)).*Xp;%calculation for rectangular parallel part
                        xp1=vpa(xp1);
                        xp2=@(x)(((u.*(h+(tand(alpha))*x))./(2.*m_bit2*w)).*Xp);%calculation for rectangular flared part
                        xp2=vpa(int(xp2,x,[0,flare_l]));
                        xp=xp1+xp2;
                        ind_grad1=(0.6+0.4.*log(h./(w+d)))*1e-6;%calculation for rectangular parallel part
                        ind_grad2=@(x)(0.6+0.4.*log((h+tand(alpha)*x)./(w+d)))*1e-6;
                        ind_grad2=vpa(int(ind_grad2,x,[0,flare_l]));%calculation for rectangular flared part
                        ind_grad=ind_grad1+ind_grad2;
                    case 3
                        m_bit1=(sigma*d_sheath*d./(lambda.*h))*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for rectangular parallel part
                        m_bit1=vpa(m_bit1);
                        m_bit2=@(x)((sigma*d_sheath.*d)./(lambda.*(h+(tand(alpha))*x))).*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for flared part
                        m_bit2=vpa(int(m_bit2,x,[0,flare_l]));
                        mbit=m_bit1+m_bit2;
                        xp1=(u.*h./(2.*m_bit1.*w)).*Xp;%calculation for rectangular parallel part
                        flare_l_a=0.8.*flare_l;%part with electrode in rectangular shape
                        flare_l_b=0.2.*flare_l;%part with electrode in tongue shape
                        xp2=@(x)(((u.*(h+(tand(alpha))*x))./(2.*m_bit2.*w)).*Xp);
                        xp2=vpa(int(xp2,x,[0,flare_l_a]));
                        xp3=@(x)(((u.*(h+(tand(alpha))*x))./(2.*m_bit2.*(w*(1-x/l)))).*Xp);
                        xp3=vpa(int(xp3,x,[0,flare_l_b]));
                        xp=xp1+xp2+xp3;
                        ind_grad1=(0.6+0.4.*log(h./(w+d)))*1e-6;%calculation for rectangular parallel part
                        ind_grad2=@(x)(0.6+0.4.*log((h+tand(alpha)*x)./(w+d)))*1e-6;%calculation for rectangular flared part
                        ind_grad2=vpa(int(ind_grad2,x,[0,flare_l_a]));
                        ind_grad3=@(x)(0.6+0.4.*log((h+tand(alpha)*x)./((w.*(1-x/l))+d)))*1e-6;%calculation for tongue flared part
                        ind_grad3=vpa(int(ind_grad3,x,[0,flare_l_b]));
                        ind_grad=ind_grad1+ind_grad2+ind_grad3;
                end                 
            end
            ve=xp;%exhaust velocity
            Isp = ve/g;
            Y=1.3;%for teflon
            impulse_bit_EM=Isp*mbit*g*f; %impulse bit
            impulse_bit_ET=sqrt((8*(Y-1)*mbit*E)/(Y^2.*(Y+1)))+(Xp*ind_grad)/2;%impulse bit
            impulse_bit= impulse_bit_EM+ impulse_bit_ET;
            Efficiency=impulse_bit.*ve/(2*E);%efficiency
            thrust=impulse_bit*f;
            w=w.*1e3;
            switch get(handles.dependent_variable_popup_large,'Value')
                case 2
                    plot(w,Efficiency*1e2,'-*');
                    xlabel('electrode width(mm)');
                    ylabel('Thruster Efficiency(%)');
                case 3
                    plot(w,impulse_bit*1e6,'-*');
                    xlabel('electrode width(mm)');
                    ylabel('Impulse Bit(uNs)');
                case 4
                    plot(w,ve,'-*');
                    xlabel('electrode width(mm)');
                    ylabel('Exhaust Velocity(m/s)');
                case 5
                    plot(w,thrust*1e6,'-*');
                    xlabel('electrode width(mm)');
                    ylabel('Thrust(uN)');
            end
    end
end
waitbar(1,progress_f,'Finished');
grid on;
pause(0.5);
close(progress_f);
set(handles.plot_button_large,'BackgroundColor',[0.85,0.85,0.85]);
guidata(hObject,handles);

function dependent_variable_popup_large_Callback(hObject, eventdata, handles)

function dependent_variable_popup_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function independent_variable_popup_large_Callback(hObject, eventdata, handles)

function independent_variable_popup_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function range_min_large_Callback(hObject, eventdata, handles)

function range_min_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function range_max_large_Callback(hObject, eventdata, handles)

function range_max_large_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function resolution_large_Callback(hObject, eventdata, handles)

function resolution_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_popup_large_Callback(hObject, eventdata, handles)

function plot_popup_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function E_A_edit_large_Callback(hObject, eventdata, handles)

function E_A_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function aspect_ratio_edit_large_Callback(hObject, eventdata, handles)

function aspect_ratio_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function exposed_propellant_edit_large_Callback(hObject, eventdata, handles)

function exposed_propellant_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function capacitance_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.capacitance_slider_large,'value',str2num(edit));
cap = str2num(get(handles.capacitance_edit_large,'String'));
minSliderValue = get(handles.capacitance_slider_large, 'Min');
maxSliderValue = get(handles.capacitance_slider_large, 'Max');
if cap < minSliderValue
    c = minSliderValue;
    set(handles.capacitance_edit_large,'String',num2str(c));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.capacitance_slider_large, 'Value', c);
elseif  cap > maxSliderValue 
    c = maxSliderValue;
    set(handles.capacitance_edit_large,'String',num2str(c));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.capacitance_slider_large, 'Value', c);
end
guidata(hObject,handles);

function capacitance_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_inductance_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.electrode_inductance_slider_large,'value',str2num(edit));
res = str2num(get(handles.electrode_inductance_edit_large,'String'));
minSliderValue = get(handles.electrode_inductance_slider_large, 'Min');
maxSliderValue = get(handles.electrode_inductance_slider_large, 'Max');
if res < minSliderValue
     r = minSliderValue ;
    set(handles.electrode_inductance_edit_large,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_inductance_slider_large, 'Value', r);
elseif  res > maxSliderValue
    r = maxSliderValue;
    set(handles.electrode_inductance_edit_large,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_inductance_slider_large, 'Value', r);
end
guidata(hObject,handles);

function electrode_inductance_edit_large_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_inductance_slider_large_Callback(hObject, eventdata, handles)
resistance=get(hObject,'value');
set(handles.electrode_inductance_edit_large,'string',num2str(resistance));
guidata(hObject,handles);

function electrode_inductance_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function capacitor_inductance_slider_large_Callback(hObject, eventdata, handles)
esr=get(hObject,'value');
set(handles.capacitor_inductance_edit_large,'string',num2str(esr));
guidata(hObject,handles);

function capacitor_inductance_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function capacitor_inductance_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.capacitor_inductance_slider_large,'value',str2num(edit));
res = str2num(get(handles.capacitor_inductance_edit_large,'String'));
minSliderValue = get(handles.capacitor_inductance_slider_large, 'Min');
maxSliderValue = get(handles.capacitor_inductance_slider_large, 'Max');
if res < minSliderValue
    r = minSliderValue ;
    set(handles.capacitor_inductance_edit_large,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.capacitor_inductance_slider_large, 'Value', r);
elseif  res > maxSliderValue 
    r = maxSliderValue;
    set(handles.capacitor_inductance_edit_large,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.capacitor_inductance_slider_large, 'Value', r);
end
guidata(hObject,handles);

function capacitor_inductance_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_resistance_slider_large_Callback(hObject, eventdata, handles)
ztot=get(hObject,'value');
set(handles.electrode_resistance_edit_large,'string',num2str(ztot));
guidata(hObject,handles);

function electrode_resistance_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function electrode_resistance_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.electrode_resistance_slider_large,'value',str2num(edit));
res = str2num(get(handles.electrode_resistance_edit_large,'String'));
minSliderValue = get(handles.electrode_resistance_slider_large, 'Min');
maxSliderValue = get(handles.electrode_resistance_slider_large, 'Max');
if res < minSliderValue
    r = minSliderValue ;
    set(handles.electrode_resistance_edit_large,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_resistance_slider_large, 'Value', r);
elseif  res > maxSliderValue
    r = maxSliderValue;
    set(handles.electrode_resistance_edit_large,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_resistance_slider_large, 'Value', r);
end
guidata(hObject,handles);

function electrode_resistance_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function capacitor_resistance_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.capacitor_resiatnce_slider_large,'value',str2num(edit));
len = str2num(get(handles.capacitor_resistance_edit_large,'String'));
minSliderValue = get(handles.capacitor_resiatnce_slider_large, 'Min');
maxSliderValue = get(handles.capacitor_resiatnce_slider_large, 'Max');
if len < minSliderValue
     l = minSliderValue ;
    set(handles.capacitor_resistance_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.capacitor_resiatnce_slider_large, 'Value', l);
elseif  len > maxSliderValue 
    l = maxSliderValue;
    set(handles.capacitor_resistance_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.capacitor_resiatnce_slider_large, 'Value', l);
end
guidata(hObject,handles);

function capacitor_resistance_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function capacitor_resiatnce_slider_large_Callback(hObject, eventdata, handles)
i=get(hObject,'value');
set(handles.capacitor_resistance_edit_large,'string',num2str(i));
guidata(hObject,handles);

function capacitor_resiatnce_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function propellant_feed_popup_large_Callback(hObject, eventdata, handles)

function propellant_feed_popup_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function propellant_width_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.propellant_width_slider_large,'value',str2num(edit));
len = str2num(get(handles.propellant_width_edit_large,'String'));
minSliderValue = get(handles.propellant_width_slider_large, 'Min');
maxSliderValue = get(handles.propellant_width_slider_large, 'Max');
if len < minSliderValue
    l = minSliderValue;
    set(handles.propellant_width_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.propellant_width_slider_large, 'Value', l);
elseif  len > maxSliderValue 
    l = maxSliderValue;
    set(handles.propellant_width_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.propellant_width_slider_large, 'Value', l);
end
guidata(hObject,handles);

function propellant_width_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function propellant_width_slider_large_Callback(hObject, eventdata, handles)
length=get(hObject,'value');
set(handles.propellant_width_edit_large,'string',num2str(length));
guidata(hObject,handles);

function propellant_width_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function propellant_length_slider_large_Callback(hObject, eventdata, handles)
length=get(hObject,'value');
set(handles.propellant_length_edit_large,'string',num2str(length));
guidata(hObject,handles);

function propellant_length_slider_large_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function propellant_length_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.propellant_length_slider_large,'value',str2num(edit));
len = str2num(get(handles.propellant_length_edit_large,'String'));
minSliderValue = get(handles.propellant_length_slider_large, 'Min');
maxSliderValue = get(handles.propellant_length_slider_large, 'Max');
if len < minSliderValue
    l = minSliderValue ;
    set(handles.propellant_length_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.propellant_length_slider_large, 'Value', l);
elseif  len > maxSliderValue 
    l = maxSliderValue;
    set(handles.propellant_length_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.propellant_length_slider_large, 'Value', l);
end
guidata(hObject,handles);

function propellant_length_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function propellant_height_slider_large_Callback(hObject, eventdata, handles)
length=get(hObject,'value');
set(handles.propellant_height_edit_large,'string',num2str(length));
guidata(hObject,handles);

function propellant_height_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function propellant_height_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.propellant_height_slider_large,'value',str2num(edit));
len = str2num(get(handles.propellant_height_edit_large,'String'));
minSliderValue = get(handles.propellant_height_slider_large, 'Min');
maxSliderValue = get(handles.propellant_height_slider_large, 'Max');
if len < minSliderValue
    l = minSliderValue ;
    set(handles.propellant_height_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.propellant_height_slider_large, 'Value', l);
elseif  len > maxSliderValue 
    l = maxSliderValue;
    set(handles.propellant_height_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.propellant_height_slider_large, 'Value', l);
end
guidata(hObject,handles);

function propellant_height_edit_large_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function flare_button_large_Callback(hObject, eventdata, handles)
handlesArray = [handles.flare_angle_slider_large,handles.flare_angle_edit_large,
    handles.flare_length_slider_large,handles.flare_length_edit_large];
if get(hObject,'Value') == 1
    set(handlesArray, 'Enable', 'on');
else
    set(handlesArray, 'Enable', 'off');
end

function calculate_final_large_Callback(hObject, eventdata, handles)
set(handles.calculate_final_large,'BackgroundColor',[1,1,1]);
digits(7);
progress_f= waitbar(0,'Please wait...');
waitbar(.33,progress_f,'Loading your data');
pause(0.5);
%calculations
u=1.2566370614*(1e-6);%permeability of free space
g=9.8;
E=str2num(get(handles.discharge_energy_edit_large,'string'));%discharge energy in J
f=str2num(get(handles.pulse_frequency_edit_large,'string'));%pulse frequency in Hz
C=(str2num(get(handles.capacitance_edit_large,'string')))*1e-6;%capacitance in F
Le=(str2num(get(handles.electrode_inductance_edit_large,'string')))*1e-9;%electrode inductance in H
Re=(str2num(get(handles.electrode_resistance_edit_large,'string')))*1e-3;%electrode resistance in ohm                              
ESR=(str2num(get(handles.capacitor_resistance_edit_large,'string')))*1e-3;%capacitor resistance in ohm 
Lcap=(str2num(get(handles.capacitor_inductance_edit_large,'string')))*1e-9;%capacitor inductance in H
h=(str2num(get(handles.electrode_spacing_edit_large,'string')))*1e-3;%electrode spacing
w=(str2num(get(handles.electrode_width_edit_large,'string')))*1e-3;%electrode width
l=(str2num(get(handles.electrode_length_edit_large,'string')))*1e-3;%electrode length
d=(str2num(get(handles.electrode_thickness_edit_large,'string')))*1e-3; %electrode thickness
t_span=(str2num(get(handles.time_edit_large,'string')))*1e-6; %time_span for plasma sheet to reach the end of electrodes                                
if (h/w <2)
    warndlg('It seems the aspect ratio is less than 1. This makes your calculation into the catrgory of Micro PPT. Kindly, check the dimensions again or use the Micro PPT caculation page.');
else  
    L_large= 0.4*l*(log(h/(w+d)) + 3/2 - h/l +0.22*(w+d)/l);
    L=(L_large*1e-6)+Lcap+Le;%total initial resistance
    R=ESR+Re;%total initial resistance
    Vo=sqrt(2*E/C);%voltage
    zi=(R/2)*sqrt(C/L);%zeta
    wo=1/sqrt(L*C);
    syms t;%declaring t as a variaable
    if zi<1
        i2=@(t)((Vo^2.*exp(-2*zi*wo*t)*(sin(sqrt(1-zi^2)*wo*t))^2)/(L^2.*(1-zi^2).*wo^2));%current equation
        i=@(t)((Vo.*exp(-zi*wo*t)*(sin(sqrt(1-(zi^2)).*wo*t)))/(L.*sqrt(1-(zi^2)).*wo));
        v=@(t)((Vo.*zi.*exp(-zi*wo*t).*sin(sqrt(1-(zi^2)).*wo*t)./sqrt(1-(zi^2)))+(Vo.*exp(-zi*wo*t).*cos(sqrt(1-(zi^2)).*wo*t))); %voltage equation
    else
        if C==(4*L/(R^2))
            i2=@(t)(Vo^2/L^2)*(t^2.*exp(-2*wo*t));%current equation
            i=@(t)(Vo/L)*(t.*exp(-wo*t));
            v=@(t)(Vo.*exp(-wo*t).*(1+wo*t));%voltage equation
        else
            i2=@(t)(Vo^2.*(exp(-2*zi*wo*t))*(sinh(wo*t*sqrt(zi^2-1)))^2)/(L^2.*(zi^2-1).*wo^2);%current equation
            i=@(t)(Vo.*(exp(-zi*wo*t))*(sinh(wo*t*sqrt((zi^2)-1))))/(L.*sqrt((zi^2)-1).*wo);
            v=@(t)((Vo.*zi.*exp(-zi*wo*t).*sinh(sqrt((zi^2)-1).*wo*t)./sqrt((zi^2)-1))+(Vo.*exp(-zi*wo*t).*cosh(sqrt((zi^2)-1).*wo*t)));
        end
    end
    V=v(t_span);
    I=i(t_span);%loop current flowing through the circuit
    I_arc=Vo*sqrt(C/L);%arc current flowing through the circuit loop
    i2_integration=int(i2,t,[0,t_span]);%evaluation of integral of (i(t))^2
    Xp=vpa(i2_integration);%making the above answer definite
    %calculation of mbit
    me=9.1.*1e-31;%mass of electron
    e=1.6.*1e-19;%elementary charge
    eo=8.85.*1e-12;%Vacuum permittivity of free space(F/m)
    Ro=(1/3)*sqrt(C/L);
    Lp=0;%plasma inductance
    Te=11605;%temperature(in K)(1 eV)
    r=(str2num(get(handles.spot_radius_large,'string')))*1e-3;%radius of the cathode spot
    S_spot=pi.*r^2;%Surface area of the cathode spot near the mixing region
    I_spot=(str2num(get(handles.spot_current_large,'string')));%For copper electrodes it has been experimentally observed that the current per observed spot
    So=S_spot*I_arc/I_spot;%Surface area of the initial plasma flow near the mixing region
    del_T= 8*11605;%temperature difference(8 eV)(in K)
    lambda=1.2*del_T;%latent heat of vaporization, specific heat=1.2
    alp_I=(str2num(get(handles.alpha_large,'string')));%Ion current normalised by arc current, the ratio between the ion current and the arc current in a cathodic plasma
    Q=(str2num(get(handles.mean_ion_charge_large,'string')));%mean ion charge state number in plasma flow
    V_spot=(str2num(get(handles.spot_velocity_large,'string')))*1e3;%initial copper plasma has a velocity towards the direction of the anode
    Ne=I*alp_I./(Q*e*So*V_spot);%Electron density
    Ni=Ne/Q;%ion density
    A=23-log(Ne^0.5*Q*Te^-1.5);%Coulomb logarithm, factor by which small angle collisions are more effective than large angle collision
    neu=(3.62.*1e-6)*A*Ne*(Te^(-3/2))/Q;%The frequency at which electrons and ions collide
    sigma=(Ne.*(e^2))/(me*neu);%conductivity
    Rp=abs(V/I);%plasma resistance
    V_sheath=abs(V/2);
    d_sheath=sqrt(2*V_sheath*eo/(e*Ni));%thickness of the sheet
    I_sheath=pi*(Re^2)*(2.33e-6).*((V_sheath^1.5)/(d_sheath^2));
    R_sheath=V_sheath/I_sheath;%sheath resistance
    Ztot=R+Rp+R_sheath;%total resistance of the circuit
    waitbar(.67,progress_f,'Processing your data');
    pause(0.5);
    syms x;
    %for non flared geometry
    if (get(handles.flare_button_large,'value'))==0
        m_bit=(sigma*d_sheath*d/(lambda*h))*(Rp^2+Lp^2/L*C)*(E/Ro);%m_ablated = m_bit
        mbit=vpa(m_bit);
        switch get(handles.electrode_shape_popup_large,'Value')
            case 2
                xp=(u*h/(2.*mbit.*w)).*Xp;
                ind_grad=(0.6+0.4*log(h/(w+d)))*1e-6;
            case 3
                tl=0.2*l;%part with electrode in tongue shape
                xp1=(u*h/(2.*mbit.*w)).*Xp;%calculation for rectangular part
                xp1=vpa(xp1);
                xp2=@(x)((u*h/(2*mbit*w.*(1-x/l)))*Xp);%calculation for tongue part
                xp2=vpa(int(xp2,x,[0,tl]));
                xp=xp1+xp2;
                ind_grad=(0.6+0.4*log(h/(w+d)))*1e-6;
        end
        ve=xp;%exhaust velocity
        Isp = ve/g;%Specific Impulse
        Y=1.3;%for teflon
        impulse_bit_EM=Isp*mbit*g*f; %impulse bit
        impulse_bit_ET=sqrt((8*(Y-1)*mbit*E)/(Y^2.*(Y+1)))+(Xp*ind_grad)/2;%impulse bit
        impulse_bit=impulse_bit_EM+impulse_bit_ET;
        Efficiency=impulse_bit.*ve/(2*E);%efficiency
        thrust=impulse_bit*f;%thrust
        tr_eff=1-(ESR/Ztot);%transfer efficiency
        acc_eff=Efficiency/tr_eff;%acceleration efficiency
        hprop=(str2num(get(handles.propellant_height_edit_large,'string')))*1e-3;
        wprop=(str2num(get(handles.propellant_width_edit_large,'string')))*1e-3;
        aspect_ratio=h/w;
        set(handles.aspect_ratio_large,'string',num2str(aspect_ratio));
        switch get(handles.propellant_feed_popup_large,'Value')
            case 2
                n=1;
            case 3
                n=2;
            case 4
                n=3;
        end
        expo_area=(n*hprop*wprop)*1e4;
        E_A=E./expo_area;
        gamma=(str2num(get(handles.ion_erosion_large,'string')));%ion erosion rate(in µg/C)
        %considering the number of pulses to be 5e5
        dm_dt=I.*gamma/5e5;%rate of mass loss from the electrode surface as a function of the discharge current
        set(handles.electrode_erosion_large,'string',num2str(dm_dt));
        set(handles.exposed_propellant_large,'string',num2str(expo_area));
        set(handles.E_A_ratio_large,'string',num2str(E_A));
        set(handles.specific_impulse_large,'string',char(Isp));
        set(handles.exhaust_velocity_large,'string',char(ve));
        set(handles.inductance_gradient_large,'string',num2str(ind_grad*1e6));
        set(handles.impulse_bit_large,'string',char(impulse_bit*1e6));
        set(handles.thruster_eff_large,'string',char(Efficiency*1e2));
        set(handles.thrust_large,'string',char(thrust*1e6));
        set(handles.transfer_efficiency_large,'string',num2str(tr_eff*1e2));
        set(handles.acceleration_eff_large,'string',char(acc_eff*1e2));
        set(handles.inductance_change_large,'string','nan');
    end
    %for flared geometry
    if (get(handles.flare_button_large,'value'))==1
        flare_l= ((str2num(get(handles.flare_length_edit_large,'string')))*1e-3);
        alpha=str2num(get(handles.flare_angle_edit_large,'string'));%flare angle
        switch get(handles.electrode_shape_popup_large,'Value')
            case 2
                m_bit1=(sigma*d_sheath*d/(lambda*h))*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for rectangular parallel part
                m_bit1=vpa(m_bit1);
                m_bit2=@(x)((sigma*d_sheath*d/(lambda.*(h+(tand(alpha))*x))).*(Rp^2+Lp^2./(L*C)).*(E/Ro));%calculation for rectangular flared part
                m_bit2=vpa(int(m_bit2,x,[0,flare_l]));
                mbit=m_bit1+m_bit2;
                xp1=(u*h/(2*m_bit1*w))*Xp;%calculation for rectangular parallel part
                xp1=vpa(xp1);
                xp2=@(x)(((u.*(h+(tand(alpha))*x))/(2*m_bit2*w))*Xp);%calculation for rectangular flared part
                xp2=vpa(int(xp2,x,[0,flare_l]));
                xp=xp1+xp2;
                ind_grad1=(0.6+0.4*log(h/(w+d)))*1e-6;%calculation for rectangular parallel part
                ind_grad2=@(x)(0.6+0.4*log((h+tand(alpha)*x)/(w+d)))*1e-6;
                ind_grad2=vpa(int(ind_grad2,x,[0,flare_l]));%calculation for rectangular flared part
                ind_grad=ind_grad1+ind_grad2;
                %calculation of inductance change
                fun_x=@(x)(u/(2*pi))*((2*(h+(tand(alpha))*x)/w)*(pi-2*atan((h+tand(alpha)*x)/w)+(h+tand(alpha)*x)/w.*log((h+tand(alpha)*x)))-2*log(w)+(1-((h+tand(alpha)*x)/w).^2)*(log((h+tand(alpha)*x).^2+w.^2)));
            case 3
                m_bit1=(sigma*d_sheath*d/(lambda*h))*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for rectangular parallel part
                m_bit1=vpa(m_bit1);
                m_bit2=@(x)((sigma*d_sheath.*d)./(lambda.*(h+(tand(alpha))*x)))*(Rp^2+Lp^2/L*C)*(E/Ro);%calculation for flared part
                m_bit2=vpa(int(m_bit2,x,[0,flare_l]));
                mbit=m_bit1+m_bit2;
                xp1=(u*h/(2*m_bit1.*w))*Xp;%calculation for rectangular parallel part
                flare_l_a=0.8.*flare_l;%part with electrode in rectangular shape
                flare_l_b=0.2.*flare_l;%part with electrode in tongue shape
                xp2=@(x)(((u.*(h+(tand(alpha))*x))/(2*m_bit2.*w))*Xp);
                xp2=vpa(int(xp2,x,[0,flare_l_a]));
                xp3=@(x)(((u.*(h+(tand(alpha))*x))/(2*m_bit2.*(w*(1-x/l))))*Xp);
                xp3=vpa(int(xp3,x,[0,flare_l_b]));
                xp=xp1+xp2+xp3;
                ind_grad1=(0.6+0.4*log(h/(w+d)))*1e-6;%calculation for rectangular parallel part
                ind_grad2=@(x)(0.6+0.4*log((h+tand(alpha)*x)/(w+d)))*1e-6;%calculation for rectangular flared part
                ind_grad2=vpa(int(ind_grad2,x,[0,flare_l_a]));
                ind_grad3=@(x)(0.6+0.4*log((h+tand(alpha)*x)/((w*(1-x/l))+d)))*1e-6;%calculation for tongue flared part
                ind_grad3=vpa(int(ind_grad3,x,[0,flare_l_b]));
                ind_grad=ind_grad1+ind_grad2+ind_grad3;
                %calculation of inductance change
                fun_x=@(x)((u/(2*pi))*((2*(h+tand(alpha)*x)/(w*(1-x/l)))*(pi-2*atan((h+tand(alpha)*x)/(w*(1-x/l)))+(h+tand(alpha)*x)/(w*(1-x/l)).*log((h+tand(alpha)*x)))-2*log((w*(1-x/l)))+(1-((h+tand(alpha)*x)/(w*(1-x/l))).^2)*(log((h+tand(alpha)*x).^2+(w*(1-x/l)).^2))));
        end
        del_L=int(fun_x,x,[0,l]);
        delL=vpa(del_L,5);
        ve=xp;%exhaust velocity
        Isp = ve/g;%Specific Impulse
        Y=1.3; %for teflon
        impulse_bit_EM=Isp*mbit*g*f; %impulse bit electromagnetic
        impulse_bit_ET=sqrt((8*(Y-1)*mbit*E)/(Y^2.*(Y+1)))+(Xp*ind_grad)/2; %impulse bit electrothermal
        impulse_bit=impulse_bit_EM+impulse_bit_ET;
        Efficiency=impulse_bit.*ve/(2*E);%efficiency
        thrust=impulse_bit*f;%thrust
        tr_eff=1-(ESR/Ztot);%transfer efficiency
        acc_eff=Efficiency/tr_eff;%acceleration efficiency
        hprop=(str2num(get(handles.propellant_height_edit_large,'string')))*1e-3;
        wprop=(str2num(get(handles.propellant_width_edit_large,'string')))*1e-3;
        aspect_ratio=h/w;
        set(handles.aspect_ratio_large,'string',num2str(aspect_ratio));
        switch get(handles.propellant_feed_popup_large,'Value')
            case 2
                n=1;
            case 3
                n=2;
            case 4
                n=3;
        end
        expo_area=(n*hprop*wprop)*(1e4);
        E_A=E./expo_area;
        gamma=(str2num(get(handles.ion_erosion_large,'string')));%ion erosion rate
        %considering the number of pulses to be 5e5
        dm_dt=I.*gamma/5e5;%rate of mass loss from the electrode surface as a function of the discharge current
        set(handles.electrode_erosion_large,'string',num2str(dm_dt));
        set(handles.exposed_propellant_large,'string',num2str(expo_area));
        set(handles.E_A_ratio_large,'string',num2str(E_A));
        set(handles.specific_impulse_large,'string',char(Isp));
        set(handles.exhaust_velocity_large,'string',char(ve));
        set(handles.inductance_gradient_large,'string',char(ind_grad*1e6));
        set(handles.impulse_bit_large,'string',char(impulse_bit*1e6));
        set(handles.thruster_eff_large,'string',char(Efficiency*1e2));
        set(handles.thrust_large,'string',char(thrust*1e6));
        set(handles.transfer_efficiency_large,'string',num2str(tr_eff*1e2));
        set(handles.acceleration_eff_large,'string',char(acc_eff*1e2));
        set(handles.inductance_change_large,'string',char(delL*1e9));
    end
end
waitbar(1,progress_f,'Finished');
pause(0.5);
close(progress_f);
set(handles.calculate_final_large,'BackgroundColor',[0.85,0.85,0.85]);
guidata(hObject,handles);

function time_edit_large_Callback(hObject, eventdata, handles)

function time_edit_large_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function flare_angle_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.flare_angle_slider_large,'value',str2num(edit));
angle = str2num(get(handles.flare_angle_edit_large,'String'));
minSliderValue = get(handles.flare_angle_slider_large, 'Min');
maxSliderValue = get(handles.flare_angle_slider_large, 'Max');
if angle < minSliderValue
    l = minSliderValue ;
    set(handles.flare_angle_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.flare_angle_slider_large, 'Value', l);
elseif  angle > maxSliderValue 
    l = maxSliderValue;
    set(handles.flare_angle_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.flare_angle_slider_large, 'Value', l);
end
guidata(hObject,handles);

function flare_angle_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function flare_angle_slider_large_Callback(hObject, eventdata, handles)
fangle=get(hObject,'value');
set(handles.flare_angle_edit_large,'string',num2str(fangle));
guidata(hObject,handles);

function flare_angle_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function flare_length_slider_large_Callback(hObject, eventdata, handles)
width=get(hObject,'value');
set(handles.flare_length_edit_large,'string',num2str(width));
guidata(hObject,handles);

function flare_length_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function flare_length_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.flare_length_slider_large,'value',str2num(edit));
angle = str2num(get(handles.flare_length_edit_large,'String'));
minSliderValue = get(handles.flare_length_slider_large, 'Min');
maxSliderValue = get(handles.flare_length_slider_large, 'Max');
if angle < minSliderValue
    l = minSliderValue ;
    set(handles.flare_length_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.flare_length_slider_large, 'Value', l);
elseif  angle > maxSliderValue 
    l = maxSliderValue;
    set(handles.flare_length_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.flare_length_slider_large, 'Value', l);
end
guidata(hObject,handles);

function flare_length_edit_large_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_spacing_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.electrode_spacing_slider_large,'value',str2num(edit));
x = str2num(get(handles.electrode_spacing_edit_large,'String'));
minSliderValue = get(handles.electrode_spacing_slider_large, 'Min');
maxSliderValue = get(handles.electrode_spacing_slider_large, 'Max');
if x < minSliderValue
    s = minSliderValue;
    set(handles.electrode_spacing_edit_large,'String',num2str(s));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_spacing_slider_large, 'Value', s);
elseif  x > maxSliderValue 
     s = maxSliderValue;
    set(handles.electrode_spacing_edit_large,'String',num2str(s));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_spacing_slider_large, 'Value', s);
end
guidata(hObject,handles);

function electrode_spacing_edit_large_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_width_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.electrode_width_slider_large,'value',str2num(edit));
x = str2num(get(handles.electrode_width_edit_large,'String'));
minSliderValue = get(handles.electrode_width_slider_large, 'Min');
maxSliderValue = get(handles.electrode_width_slider_large, 'Max');
if x < minSliderValue
     r = minSliderValue ;
    set(handles.electrode_width_edit_large,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_width_slider_large, 'Value', r);
elseif  x > maxSliderValue 
    r = maxSliderValue;
    set(handles.electrode_width_edit_large,'String',num2str(r));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_width_slider_large, 'Value', r);
end
guidata(hObject,handles);

function electrode_width_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_width_slider_large_Callback(hObject, eventdata, handles)
width=get(hObject,'value');
set(handles.electrode_width_edit_large,'string',num2str(width));
guidata(hObject,handles);

function electrode_width_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function electrode_spacing_slider_large_Callback(hObject, eventdata, handles)
spacing=get(hObject,'value');
set(handles.electrode_spacing_edit_large,'string',num2str(spacing));
guidata(hObject,handles);

function electrode_spacing_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function electrode_length_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.electrode_length_slider_large,'value',str2num(edit));
len = str2num(get(handles.electrode_length_edit_large,'String'));
minSliderValue = get(handles.electrode_length_slider_large, 'Min');
maxSliderValue = get(handles.electrode_length_slider_large, 'Max');
if len < minSliderValue
     l = minSliderValue;
    set(handles.electrode_length_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
       set(handles.electrode_length_slider_large, 'Value', l);
elseif  len > maxSliderValue 
     l = maxSliderValue;
    set(handles.electrode_length_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
       set(handles.electrode_length_slider_large, 'Value', l);
end
guidata(hObject,handles);

function electrode_length_edit_large_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_length_slider_large_Callback(hObject, eventdata, handles)
length=get(hObject,'value');
set(handles.electrode_length_edit_large,'string',num2str(length));
guidata(hObject,handles);

function electrode_length_slider_large_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function electrode_shape_popup_large_Callback(hObject, eventdata, handles)

function electrode_shape_popup_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_thickness_slider_large_Callback(hObject, eventdata, handles)
thickness=get(hObject,'value');
set(handles.electrode_thickness_edit_large,'string',num2str(thickness));
guidata(hObject,handles);

function electrode_thickness_slider_large_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function electrode_thickness_edit_large_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.electrode_thickness_slider_large,'value',str2num(edit));
len = str2num(get(handles.electrode_thickness_edit_large,'String'));
minSliderValue = get(handles.electrode_thickness_slider_large, 'Min');
maxSliderValue = get(handles.electrode_thickness_slider_large, 'Max');
if len < minSliderValue
    l = minSliderValue ;
    set(handles.electrode_thickness_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_thickness_slider_large, 'Value', l);
elseif  len > maxSliderValue 
    l = maxSliderValue;
    set(handles.electrode_thickness_edit_large,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_thickness_slider_large, 'Value', l);
end
guidata(hObject,handles);

function electrode_thickness_edit_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_material_large_Callback(hObject, eventdata, handles)

function electrode_material_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_material_Callback(hObject, eventdata, handles)

function electrode_material_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ion_erosion_Callback(hObject, eventdata, handles)

function ion_erosion_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function current_per_spot_Callback(hObject, eventdata, handles)

function current_per_spot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function spot_radius_micro_Callback(hObject, eventdata, handles)

function spot_radius_micro_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function spot_velocity_Callback(hObject, eventdata, handles)

function spot_velocity_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Mean_ion_charge_Callback(hObject, eventdata, handles)

function Mean_ion_charge_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ion_erosion_large_Callback(hObject, eventdata, handles)

function ion_erosion_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function spot_current_large_Callback(hObject, eventdata, handles)

function spot_current_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function spot_radius_large_Callback(hObject, eventdata, handles)

function spot_radius_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function spot_velocity_large_Callback(hObject, eventdata, handles)

function spot_velocity_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mean_ion_charge_large_Callback(hObject, eventdata, handles)

function mean_ion_charge_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dm_dt_Callback(hObject, eventdata, handles)

function dm_dt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function alpha_i_Callback(hObject, eventdata, handles)

function alpha_i_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_thickness_sld_Callback(hObject, eventdata, handles)
thickness=get(hObject,'value');
set(handles.electrode_thickness_edit,'string',num2str(thickness));
guidata(hObject,handles);

function electrode_thickness_sld_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function electrode_thickness_edit_Callback(hObject, eventdata, handles)
edit=get(hObject,'string');
set(handles.electrode_thickness_sld,'value',str2num(edit));
len = str2num(get(handles.electrode_thickness_edit,'String'));
minSliderValue = get(handles.electrode_thickness_sld, 'Min');
maxSliderValue = get(handles.electrode_thickness_sld, 'Max');
if len < minSliderValue
    l = minSliderValue ;
    set(handles.electrode_thickness_edit,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_thickness_sld, 'Value', l);
elseif  len > maxSliderValue 
    l = maxSliderValue;
    set(handles.electrode_thickness_edit,'String',num2str(l));
    warndlg('Value Out of Range! Please Enter Again');
    set(handles.electrode_thickness_sld, 'Value', l);
end
guidata(hObject,handles);

function electrode_thickness_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function alpha_large_Callback(hObject, eventdata, handles)

function alpha_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function aspect_ratio_large_Callback(hObject, eventdata, handles)

function aspect_ratio_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function exposed_propellant_large_Callback(hObject, eventdata, handles)

function exposed_propellant_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function E_A_ratio_large_Callback(hObject, eventdata, handles)

function E_A_ratio_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function electrode_erosion_large_Callback(hObject, eventdata, handles)

function electrode_erosion_large_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%the below part has been made to avoid a conflict
function popupmenu17_Callback(hObject, eventdata, handles)

function popupmenu17_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu18_Callback(hObject, eventdata, handles)

function popupmenu18_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu19_Callback(hObject, eventdata, handles)

function popupmenu19_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit168_Callback(hObject, eventdata, handles)

function edit168_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit169_Callback(hObject, eventdata, handles)

function edit169_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit170_Callback(hObject, eventdata, handles)

function edit170_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function micro_ppt_CreateFcn(hObject, eventdata, handles)

function micro_ppt_DeleteFcn(hObject, eventdata, handles)
