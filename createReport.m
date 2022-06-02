%% createReport.m
% Function to create reports from result and setup data from the
% simulation.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function createReport()

    %% Import API Packages and Define Document Type
    import mlreportgen.*;
    import mlreportgen.dom.*;

    doctype = "pdf";
    rpt = report("Laptime Simulation Report", doctype);
    
    %% Create Title Page
    tp = TitlePage("Title","Laptime Simulation Report","Subtitle","Test","Image","balticracing_logo_transparent.png","Author","Eric Dornieden - Baltic Racing","PubDate","Report Published at now");
    add(rpt, tp);

    %% Create Table of Contents
    toc = TableOfContents;
    add(rpt, toc);

    %% Create Chapter 1 and Section 1


    close(rpt);
    rptview(rpt);

end