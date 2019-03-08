// trackingDlg.cpp : implementation file
//

#include "stdafx.h"
#include "tracking.h"
#include "trackingDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CTrackingDlg dialog

CTrackingDlg::CTrackingDlg(CWnd* pParent /*=NULL*/)
    : CSimulationDialog(CTrackingDlg::IDD, pParent)
{
    m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
    m_data.params = model::make_default_parameters();
}

void CTrackingDlg::DoDataExchange(CDataExchange* pDX)
{
    CSimulationDialog::DoDataExchange(pDX);
    DDX_Control(pDX, IDC_SLIDER1, m_frameSlider);
    DDX_Control(pDX, IDC_VIDEO, m_videoCtrl);
    DDX_Text(pDX, IDC_EDIT1, m_data.params.target_resolution_w);
    DDX_Text(pDX, IDC_EDIT2, m_data.params.target_resolution_h);
}

BEGIN_MESSAGE_MAP(CTrackingDlg, CSimulationDialog)
    ON_WM_PAINT()
    ON_WM_QUERYDRAGICON()
    ON_BN_CLICKED(IDC_BUTTON1, &CTrackingDlg::OnBnClickedButton1)
    ON_BN_CLICKED(IDC_BUTTON2, &CTrackingDlg::OnBnClickedButton2)
    ON_WM_HSCROLL()
    ON_BN_CLICKED(IDC_BUTTON3, &CTrackingDlg::OnBnClickedButton3)
END_MESSAGE_MAP()

// CTrackingDlg message handlers

BOOL CTrackingDlg::OnInitDialog()
{
    CSimulationDialog::OnInitDialog();

    // Set the icon for this dialog.  The framework does this automatically
    //  when the application's main window is not a dialog
    SetIcon(m_hIcon, TRUE);            // Set big icon
    SetIcon(m_hIcon, FALSE);        // Set small icon

    // TODO: Add extra initialization here

    m_videoCtrl.Init(m_data);
    m_videoCtrl.callback = [this] (const geom::point < int > & p, bool lbutton) {
        if (lbutton) m_data.decorator.markers.push_back(p);
        else m_data.decorator.markers.clear();
        m_videoCtrl.RedrawBuffer();
        m_videoCtrl.SwapBuffers();
        m_videoCtrl.RedrawWindow();
    };

    return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CTrackingDlg::OnPaint()
{
    if (IsIconic())
    {
        CPaintDC dc(this); // device context for painting

        SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

        // Center icon in client rectangle
        int cxIcon = GetSystemMetrics(SM_CXICON);
        int cyIcon = GetSystemMetrics(SM_CYICON);
        CRect rect;
        GetClientRect(&rect);
        int x = (rect.Width() - cxIcon + 1) / 2;
        int y = (rect.Height() - cyIcon + 1) / 2;

        // Draw the icon
        dc.DrawIcon(x, y, m_hIcon);
    }
    else
    {
        CSimulationDialog::OnPaint();
    }
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CTrackingDlg::OnQueryDragIcon()
{
    return static_cast<HCURSOR>(m_hIcon);
}


void CTrackingDlg::OnBnClickedButton1()
{
    UpdateData(TRUE);
    m_data.source.frame = m_frameSlider.GetPos();
    StartSimulationThread();
}


void CTrackingDlg::OnBnClickedButton2()
{
    StopSimulationThread();
}


void CTrackingDlg::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar)
{
    // TODO: Add your message handler code here and/or call default

    CSimulationDialog::OnHScroll(nSBCode, nPos, pScrollBar);
}


void CTrackingDlg::OnBnClickedButton3()
{
    UpdateData(TRUE);
    CFileDialog fd(TRUE);
    if (fd.DoModal() == IDOK)
    {
        std::wstring path(fd.GetPathName().GetBuffer());
        std::string asciipath(path.begin(), path.end());
        cv::Size res(m_data.params.target_resolution_w, m_data.params.target_resolution_h);
        m_data.source.video.from_file(asciipath, res);
    }
    if (!m_data.source.video.frames.empty()) m_data.source.frame = 0;
    m_frameSlider.SetPos(0);
    m_frameSlider.SetRange(0, m_data.source.video.frames.size() - 1, TRUE);
    m_videoCtrl.RedrawBuffer();
    m_videoCtrl.SwapBuffers();
    m_videoCtrl.RedrawWindow();
}

void CTrackingDlg::OnSimulation()
{
    for (; m_bWorking && (m_data.source.frame < m_data.source.video.frames.size());
         ++m_data.source.frame)
    {
        Sleep(40);
        m_videoCtrl.RedrawBuffer();
        m_videoCtrl.SwapBuffers();
        Invoke([this] () {
            m_frameSlider.SetPos(m_data.source.frame);
            m_videoCtrl.RedrawWindow();
        });
    }
}