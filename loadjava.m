disp('loading jar file...')

%projectname = 'SimJNSNew';
if exist('C:\Users\bortk','dir')
    WINDOWS_PATH=['C:\Users\bortk\Documents\NetBeansProjects\' projectname '\dist\' projectname '.jar'];
elseif exist('C:\Users\logan','dir')
    WINDOWS_PATH=['C:\Users\logan\Documents\NetBeansProjects\' projectname '\dist\' projectname '.jar'];
else
    WINDOWS_PATH=['C:\Users\Logan\Documents\NetBeansProjects\' projectname '\dist\' projectname '.jar'];
end
if exist('./dist')
    LINUX_PATH=['./dist/' projectname '.jar'];
else
    LINUX_PATH=['~/NetBeansProjects/' projectname '/dist/' projectname '.jar'];
end

if strcmp('Windows_NT',getenv('OS'))
    disp(['using windows path: ' WINDOWS_PATH]);
    javaaddpath(WINDOWS_PATH)
else
    disp(['using linux path: ' LINUX_PATH]);
    javaaddpath(GetFullPath(LINUX_PATH));
end
