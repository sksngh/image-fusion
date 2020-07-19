%PROBLEM:Multi-modality image fusion: Every single modality image has its own
%drawbacks in providing needed information because each image is captured
%with different radiation power. In order to overcome this it is highly
%required to obtain information from multiple modalities which is used for
%clinical diagnosis. In this situatio n, fusion is a technique used to
%combine multimodality medical images such as CT, MRI, PET etc.

function varargout = fusion_gui(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @fusion_gui_OpeningFcn, ...
    'gui_OutputFcn',  @fusion_gui_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before fusion_gui is made visible.
function fusion_gui_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for fusion_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = fusion_gui_OutputFcn(hObject, eventdata, handles)
global I1;
X1=imread('1MRI.png');
X1=imresize(X1,[184 184]);
I1=rgb2gray(X1);
global I2;
X2=imread('1CT.png');
X2=imresize(X2,[184 184]);
I2=rgb2gray(X2);
I11=im2double(I1);
I22=im2double(I2);

[c1,s1]=wavedec2(I11,1,'db4');
[c2,s2]=wavedec2(I22,1,'db4');
[h1,v1,d1]=detcoef2('all',c1,s1,1);
[h2,v2,d2]=detcoef2('all',c2,s2,1);
a1 =appcoef2(c1,s1,'db4',1);
a2 =appcoef2(c2,s2,'db4',1);

X=I22;
[SIZEX,SIZEY]=size(X);


m=1.0;
delta=2^m;



N=20;
A=-1/sqrt(2*pi);
for index_x=1:N;
    for index_y=1:N;
        x=index_x-(N+1)/2;
        y=index_y-(N+1)/2;
        phi_x(index_x,index_y)=A*(x/delta^2).*exp(-(x.*x+y.*y)/(2*delta^2));
        phi_y(index_x,index_y)=A*(y/delta^2).*exp(-(x.*x+y.*y)/(2*delta^2));
    end
end;
phi_x=phi_x/norm(phi_x);
phi_y=phi_y/norm(phi_y);

Gx=conv2(X,phi_x,'same');
Gy=conv2(X,phi_y,'same');

Grads=sqrt((Gx*Gx)+(Gy*Gy));

angle_array=zeros(SIZEX,SIZEY);
for i=1:SIZEX;
    for j=1:SIZEY
        if (abs(Gx(i,j))>eps*100)
            p=atan(Gy(i,j)/Gx(i,j))*180/pi;
            if (p<0)
                p=p+360;
            end;
            if (Gx(i,j)<0 & p>180)
                p=p-180;
            elseif (Gx(i,j)<0 & p<180)
                p=p+180;
            end
        else
            p=90;
        end
        angle_array(i,j)=p;
    end
end


n=1
[c4,s4]=wavedec2(Grads,n,'db4');
[h4,v4,d4]=detcoef2('all',c4,s4,1);

n=1
[c5,s5]=wavedec2(angle_array,n,'db4');
[h5,v5,d5]=detcoef2('all',c5,s5,1);

[SIZEX,SIZEY]=size(h5);
H3=h2;
for i=2:SIZEX-1
    for j=2:SIZEY-1
        if ((h5(i,j)>=(-22.5) & h5(i,j)<=22.5) | ...
                (h5(i,j)>=(180-22.5) & h5(i,j)<=(180+22.5)))     %  0/180
            if (h4(i,j)>h4(i+1,j) & h4(i,j)>h4(i-1,j))
                H3(i,j)=h2(i,j);
            elseif h4(i+1,j)>h4(i-1,j)
                H3(i,j)=h2(i+1,j);
            else
                H3(i,j)=h2(i-1,j);
            end
            
        elseif ((h5(i,j)>=(90-22.5) & h5(i,j)<=(90+22.5)) | ...
                (h5(i,j)>=(270-22.5) & h5(i,j)<=(270+22.5))) %  90/270
            if (h4(i,j)>h4(i,j+1) & h4(i,j)>h4(i,j-1))
                H3(i,j)=h2(i,j);
            elseif h4(i,j+1)>h4(i,j-1)
                H3(i,j)=h2(i,j+1);
            else
                H3(i,j)=h2(i,j-1);
            end
            
        elseif ((h5(i,j)>=(45-22.5) & h5(i,j)<=(45+22.5)) | ...
                (h5(i,j)>=(225-22.5) & h5(i,j)<=(225+22.5))) %  45/225
            if (h4(i,j)>h4(i+1,j+1) & h4(i,j)>h4(i-1,j-1))
                H3(i,j)=h2(i,j);
            elseif h4(i+1,j+1)>h4(i-1,j-1)
                H3(i,j)=h2(i+1,j+1);
            else
                H3(i,j)=h2(i-1,j-1);
            end
            
        else  %  135/215
            if (h4(i,j)>h4(i+1,j-1) & h4(i,j)>h4(i-1,j+1))
                H3(i,j)=h2(i,j);
            elseif h4(i+1,j-1)>h4(i-1,j+1)
                H3(i,j)=h2(i+1,j-1);
            else
                H3(i,j)=h2(i-1,j+1);
            end
            
        end
    end
end




[SIZEX,SIZEY]=size(v5);
V3=v2;
for i=2:SIZEX-1
    for j=2:SIZEY-1
        if ((v5(i,j)>=(-22.5) & v5(i,j)<=22.5) | ...
                (v5(i,j)>=(180-22.5) & v5(i,j)<=(180+22.5)))     %  0/180
            if (v4(i,j)>v4(i+1,j) & v4(i,j)>v4(i-1,j))
                V3(i,j)=v2(i,j);
            elseif v4(i+1,j)>v4(i-1,j)
                V3(i,j)=v2(i+1,j);
            else
                V3(i,j)=v2(i-1,j);
            end
            
        elseif ((v5(i,j)>=(90-22.5) & v5(i,j)<=(90+22.5)) | ...
                (v5(i,j)>=(270-22.5) & v5(i,j)<=(270+22.5))) %  90/270
            if (v4(i,j)>v4(i,j+1) & v4(i,j)>v4(i,j-1))
                V3(i,j)=v2(i,j);
            elseif v4(i,j+1)>v4(i,j-1)
                V3(i,j)=v2(i,j+1);
            else
                V3(i,j)=v2(i,j-1);
            end
            
        elseif ((v5(i,j)>=(45-22.5) & v5(i,j)<=(45+22.5)) | ...
                (v5(i,j)>=(225-22.5) & v5(i,j)<=(225+22.5))) %  45/225
            if (v4(i,j)>v4(i+1,j+1) & v4(i,j)>v4(i-1,j-1))
                V3(i,j)=v2(i,j);
            elseif v4(i+1,j+1)>v4(i-1,j-1)
                V3(i,j)=v2(i+1,j+1);
            else
                V3(i,j)=v2(i-1,j-1);
            end
            
        else  %  135/215
            if (v4(i,j)>v4(i+1,j-1) & v4(i,j)>v4(i-1,j+1))
                V3(i,j)=v2(i,j);
            elseif v4(i+1,j-1)>v4(i-1,j+1)
                V3(i,j)=v2(i+1,j-1);
            else
                V3(i,j)=v2(i-1,j+1);
            end
            
        end
    end
end



[SIZEX,SIZEY]=size(d5);
D3=d2;
for i=2:SIZEX-1
    for j=2:SIZEY-1
        if ((d5(i,j)>=(-22.5) & d5(i,j)<=22.5) | ...
                (d5(i,j)>=(180-22.5) & d5(i,j)<=(180+22.5)))     %  0/180
            if (d4(i,j)>d4(i+1,j) & d4(i,j)>d4(i-1,j))
                D3(i,j)=d2(i,j);
            elseif d4(i+1,j)>d4(i-1,j)
                D3(i,j)=d2(i+1,j);
            else
                D3(i,j)=d2(i-1,j);
            end
            %        D3(i,j)=max(d23(i,j),d23(i+1,j),d23(i-1,j));
        elseif ((d5(i,j)>=(90-22.5) & d5(i,j)<=(90+22.5)) | ...
                (d5(i,j)>=(270-22.5) & d5(i,j)<=(270+22.5))) %  90/270
            if (d4(i,j)>d4(i,j+1) & d4(i,j)>d4(i,j-1))
                D3(i,j)=d2(i,j);
            elseif d4(i,j+1)>d4(i,j-1)
                D3(i,j)=d2(i,j+1);
            else
                D3(i,j)=d2(i,j-1);
            end
        elseif ((d5(i,j)>=(45-22.5) & d5(i,j)<=(45+22.5)) | ...
                (d5(i,j)>=(225-22.5) & d5(i,j)<=(225+22.5))) %  45/225
            if (d4(i,j)>d4(i+1,j+1) & d4(i,j)>d4(i-1,j-1))
                D3(i,j)=d2(i,j);
            elseif d4(i+1,j+1)>d4(i-1,j-1)
                D3(i,j)=d2(i+1,j+1);
            else
                D3(i,j)=d2(i-1,j-1);
            end
        else
            if (d4(i,j)>d4(i+1,j-1) & d4(i,j)>d4(i-1,j+1))
                D3(i,j)=d2(i,j);
            elseif d4(i+1,j-1)>d4(i-1,j+1)
                D3(i,j)=d2(i+1,j-1);
            else
                D3(i,j)=d2(i-1,j+1);
            end
        end
    end
end

Y=I11;
[SIZEX,SIZEY]=size(Y);


m=1.0;
delta=2^m;



N=20;
A=-1/sqrt(2*pi);
for index_x=1:N;
    for index_y=1:N;
        x=index_x-(N+1)/2;
        y=index_y-(N+1)/2;
        phi_x(index_x,index_y)=A*(x/delta^2).*exp(-(x.*x+y.*y)/(2*delta^2));
        phi_y(index_x,index_y)=A*(y/delta^2).*exp(-(x.*x+y.*y)/(2*delta^2));
    end
end;
phi_x=phi_x/norm(phi_x);
phi_y=phi_y/norm(phi_y);

Gx=conv2(X,phi_x,'same');
Gy=conv2(X,phi_y,'same');

Grads=sqrt((Gx*Gx)+(Gy*Gy));

angle_array=zeros(SIZEX,SIZEY);
for i=1:SIZEX;
    for j=1:SIZEY
        if (abs(Gx(i,j))>eps*100)
            p=atan(Gy(i,j)/Gx(i,j))*180/pi;
            if (p<0)
                p=p+360;
            end;
            if (Gx(i,j)<0 & p>180)
                p=p-180;
            elseif (Gx(i,j)<0 & p<180)
                p=p+180;
            end
        else
            p=90;
        end
        angle_array(i,j)=p;
    end
end


n=1
[c4,s4]=wavedec2(Grads,n,'db4');
[h4,v4,d4]=detcoef2('all',c4,s4,1);

n=1
[c5,s5]=wavedec2(angle_array,n,'db4');
[h5,v5,d5]=detcoef2('all',c5,s5,1);

[SIZEX,SIZEY]=size(h5);
H4=h1;
for i=2:SIZEX-1
    for j=2:SIZEY-1
        if ((h5(i,j)>=(-22.5) & h5(i,j)<=22.5) | ...
                (h5(i,j)>=(180-22.5) & h5(i,j)<=(180+22.5)))     %  0/180
            if (h4(i,j)>h4(i+1,j) & h4(i,j)>h4(i-1,j))
                H4(i,j)=h1(i,j);
            elseif h4(i+1,j)>h4(i-1,j)
                H4(i,j)=h1(i+1,j);
            else
                H4(i,j)=h1(i-1,j);
            end
            
        elseif ((h5(i,j)>=(90-22.5) & h5(i,j)<=(90+22.5)) | ...
                (h5(i,j)>=(270-22.5) & h5(i,j)<=(270+22.5))) %  90/270
            if (h4(i,j)>h4(i,j+1) & h4(i,j)>h4(i,j-1))
                H4(i,j)=h1(i,j);
            elseif h4(i,j+1)>h4(i,j-1)
                H4(i,j)=h1(i,j+1);
            else
                H4(i,j)=h1(i,j-1);
            end
            
        elseif ((h5(i,j)>=(45-22.5) & h5(i,j)<=(45+22.5)) | ...
                (h5(i,j)>=(225-22.5) & h5(i,j)<=(225+22.5))) %  45/225
            if (h4(i,j)>h4(i+1,j+1) & h4(i,j)>h4(i-1,j-1))
                H4(i,j)=h1(i,j);
            elseif h4(i+1,j+1)>h4(i-1,j-1)
                H4(i,j)=h1(i+1,j+1);
            else
                H4(i,j)=h1(i-1,j-1);
            end
            
        else  %  135/215
            if (h4(i,j)>h4(i+1,j-1) & h4(i,j)>h4(i-1,j+1))
                H4(i,j)=h1(i,j);
            elseif h4(i+1,j-1)>h4(i-1,j+1)
                H4(i,j)=h1(i+1,j-1);
            else
                H4(i,j)=h1(i-1,j+1);
            end
            
        end
    end
end




[SIZEX,SIZEY]=size(v5);
V4=v1;
for i=2:SIZEX-1
    for j=2:SIZEY-1
        if ((v5(i,j)>=(-22.5) & v5(i,j)<=22.5) | ...
                (v5(i,j)>=(180-22.5) & v5(i,j)<=(180+22.5)))     %  0/180
            if (v4(i,j)>v4(i+1,j) & v4(i,j)>v4(i-1,j))
                V4(i,j)=v1(i,j);
            elseif v4(i+1,j)>v4(i-1,j)
                V4(i,j)=v1(i+1,j);
            else
                V4(i,j)=v1(i-1,j);
            end
            
        elseif ((v5(i,j)>=(90-22.5) & v5(i,j)<=(90+22.5)) | ...
                (v5(i,j)>=(270-22.5) & v5(i,j)<=(270+22.5))) %  90/270
            if (v4(i,j)>v4(i,j+1) & v4(i,j)>v4(i,j-1))
                V4(i,j)=v1(i,j);
            elseif v4(i,j+1)>v4(i,j-1)
                V4(i,j)=v1(i,j+1);
            else
                V4(i,j)=v1(i,j-1);
            end
            
        elseif ((v5(i,j)>=(45-22.5) & v5(i,j)<=(45+22.5)) | ...
                (v5(i,j)>=(225-22.5) & v5(i,j)<=(225+22.5))) %  45/225
            if (v4(i,j)>v4(i+1,j+1) & v4(i,j)>v4(i-1,j-1))
                V4(i,j)=v1(i,j);
            elseif v4(i+1,j+1)>v4(i-1,j-1)
                V4(i,j)=v1(i+1,j+1);
            else
                V4(i,j)=v1(i-1,j-1);
            end
            
        else
            if (v4(i,j)>v4(i+1,j-1) & v4(i,j)>v4(i-1,j+1))
                V4(i,j)=v1(i,j);
            elseif v4(i+1,j-1)>v4(i-1,j+1)
                V4(i,j)=v1(i+1,j-1);
            else
                V4(i,j)=v1(i-1,j+1);
            end
            
        end
    end
end



[SIZEX,SIZEY]=size(d5);
D4=d1;
for i=2:SIZEX-1
    for j=2:SIZEY-1
        if ((d5(i,j)>=(-22.5) & d5(i,j)<=22.5) | ...
                (d5(i,j)>=(180-22.5) & d5(i,j)<=(180+22.5)))
            if (d4(i,j)>d4(i+1,j) & d4(i,j)>d4(i-1,j))
                D4(i,j)=d1(i,j);
            elseif d4(i+1,j)>d4(i-1,j)
                D4(i,j)=d1(i+1,j);
            else
                D4(i,j)=d1(i-1,j);
            end
        elseif ((d5(i,j)>=(90-22.5) & d5(i,j)<=(90+22.5)) | ...
                (d5(i,j)>=(270-22.5) & d5(i,j)<=(270+22.5)))
            if (d4(i,j)>d4(i,j+1) & d4(i,j)>d4(i,j-1))
                D4(i,j)=d1(i,j);
            elseif d4(i,j+1)>d4(i,j-1)
                D4(i,j)=d1(i,j+1);
            else
                D4(i,j)=d1(i,j-1);
            end
        elseif ((d5(i,j)>=(45-22.5) & d5(i,j)<=(45+22.5)) | ...
                (d5(i,j)>=(225-22.5) & d5(i,j)<=(225+22.5)))
            if (d4(i,j)>d4(i+1,j+1) & d4(i,j)>d4(i-1,j-1))
                D4(i,j)=d1(i,j);
            elseif d4(i+1,j+1)>d4(i-1,j-1)
                D4(i,j)=d1(i+1,j+1);
            else
                D4(i,j)=d1(i-1,j-1);
            end
        else
            if (d4(i,j)>d4(i+1,j-1) & d4(i,j)>d4(i-1,j+1))
                D4(i,j)=d1(i,j);
            elseif d4(i+1,j-1)>d4(i-1,j+1)
                D4(i,j)=d1(i+1,j-1);
            else
                D4(i,j)=d1(i-1,j+1);
            end
        end
    end
end

A1=(a1+a2);
H = max(H3,H4);
V = max(V3,V4);
D = max(D3,D4);
global F;
F=idwt2(A1,H,V,D,'db4');


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global I1;
axes(handles.axes1);
imshow(I1);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global I2;
axes(handles.axes2);
imshow(I2);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
global F;
axes(handles.axes3);
imshow(F);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
global I1;
global I2;
global F;
MI_fused=analysis_MI(I1,I2,F);
set(handles.text1, 'String',MI_fused );


% --- Executes when uipanel1 is resized.
function uipanel1_SizeChangedFcn(hObject, eventdata, handles)

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
global F;
global I2;
C=corr2(I2,F);
set(handles.text3, 'String',C );


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
global F;
global I1;
global I2;
C=corr2(I1,F);
set(handles.text2, 'String',C );

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
global F;
global I1;
peaksnr1 = psnr(uint8(F), uint8(I1));
set(handles.text4, 'String',peaksnr1 );

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
global F;
global I2;
peaksnr2 = psnr(uint8(F), uint8(I2));
set(handles.text5, 'String',peaksnr2);
