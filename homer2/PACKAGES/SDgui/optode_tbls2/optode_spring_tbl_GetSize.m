function tbl_size=optode_spring_tbl_GetSize(handles)

    hObject = handles.optode_spring_tbl;
    userdata = get(hObject,'userdata');
    tbl_size = userdata.tbl_size;
