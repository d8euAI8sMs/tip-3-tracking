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
    , m_bBlur(FALSE)
{
    m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
    m_data.params = model::make_default_parameters();
    m_cfg = model::tracker::make_default_config();
}

void CTrackingDlg::DoDataExchange(CDataExchange* pDX)
{
    CSimulationDialog::DoDataExchange(pDX);
    DDX_Control(pDX, IDC_SLIDER1, m_frameSlider);
    DDX_Control(pDX, IDC_VIDEO, m_videoCtrl);
    DDX_Text(pDX, IDC_EDIT1, m_data.params.target_resolution_w);
    DDX_Text(pDX, IDC_EDIT2, m_data.params.target_resolution_h);
    DDX_Text(pDX, IDC_EDIT3, m_data.params.max_frames);
    DDX_Text(pDX, IDC_EDIT4, m_cfg.mip_count);
    DDX_Text(pDX, IDC_EDIT5, m_cfg.pivot_count);
    DDX_Text(pDX, IDC_EDIT6, m_cfg.pivot_neighborhood);
    DDX_Text(pDX, IDC_EDIT7, m_cfg.eigenvalue_threshold);
    DDX_Check(pDX, IDC_CHECK1, m_bBlur);
    DDX_Check(pDX, IDC_CHECK2, m_bStrongBlur);
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
    m_videoCtrl.callback = [this] (const geom::point2d_t & p, bool lbutton) {
        if (lbutton)
        {
            if ((m_data.decorator.markers.size() & 1) == 1)
            {
                auto p0 = m_data.decorator.markers.back();
                auto p1 = p;
                if (p0.x > p1.x) std::swap(p0.x, p1.x);
                if (p0.y > p1.y) std::swap(p0.y, p1.y);
                m_data.decorator.bounding_boxes.emplace_back(
                    p0.x, p1.x, p0.y, p1.y
                );
                m_data.decorator.markers.pop_back();
            }
            else
            {
                m_data.decorator.markers.push_back(p);
            }
        }
        else
        {
            m_data.decorator.markers.clear();
            m_data.decorator.bounding_boxes.clear();
        }
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
    m_cfg.blur = (m_bBlur == TRUE);
    m_cfg.strong_blur = (m_bStrongBlur == TRUE);
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
        m_data.source.video.from_file(asciipath, m_data.params.max_frames, res);
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
    model::tracker t(m_data.params);
    for (; m_data.source.frame + 1 < m_data.source.video.frames.size();
         ++m_data.source.frame)
    {
        t.track(m_data.decorator.bounding_boxes, m_data.source, m_cfg);
        Sleep(40);
        m_videoCtrl.RedrawBuffer();
        m_videoCtrl.SwapBuffers();
        Invoke([this] () {
            m_frameSlider.SetPos(m_data.source.frame);
            m_videoCtrl.RedrawWindow();
        });
        if (!m_bWorking) break;
    }
}