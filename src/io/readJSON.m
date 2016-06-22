function info = readJSON(fname)
    fid = fopen(fname);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    info = jsondecode(str);
end