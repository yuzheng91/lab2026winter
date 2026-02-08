import React from 'react';
import { AppBar, Toolbar, Typography, Button, Box, Container } from '@mui/material';
import { useNavigate, useLocation } from 'react-router-dom'; // [新增] 引入 useLocation
import HomeIcon from '@mui/icons-material/Home';
import BuildIcon from '@mui/icons-material/Build';
import HelpIcon from '@mui/icons-material/Help';
import ContactMailIcon from '@mui/icons-material/ContactMail';

const MENU_ITEMS = [
  { label: 'Home', icon: <HomeIcon />, path: '/' },
  { label: 'Tool', icon: <BuildIcon />, path: '/wormbase' },
  { label: 'Help', icon: <HelpIcon />, path: '/help' },
  { label: 'Contact', icon: <ContactMailIcon />, path: '/contact' },
];

const Header = () => {
  const navigate = useNavigate();
  const location = useLocation(); // [步驟 1] 取得當前網址資訊

  return (
    <AppBar position="static" sx={{ bgcolor: '#222831' }}>
      <Container maxWidth="xl">
        <Toolbar disableGutters>
          <Typography
            variant="h5"
            noWrap
            component="div"
            onClick={() => navigate('/')}
            sx={{ 
              flexGrow: 1, 
              display: 'flex', 
              fontWeight: 'bold', 
              color: '#fff', 
              cursor: 'pointer' 
            }}
          >
            WGLA: Worm Gene List Analyzer
          </Typography>

          <Box sx={{ display: 'flex', gap: 1 }}>
            {MENU_ITEMS.map((item) => {
              // [步驟 2] 判斷是否為當前頁面
              // 這裡使用嚴格比對，如果你希望 /wormbase/detail 也亮燈，可以用 .startsWith
              const isActive = location.pathname === item.path;

              return (
                <Button
                  key={item.label}
                  onClick={() => navigate(item.path)}
                  sx={{
                    display: 'flex',
                    flexDirection: 'column',
                    alignItems: 'center',
                    textTransform: 'none',
                    fontSize: '0.85rem',
                    minWidth: 80,
                    py: 1,
                    transition: 'all 0.3s ease', // 讓變化有動畫效果
                    
                    // [步驟 3] 根據 isActive 設定樣式
                    // 亮燈狀態 (Active)
                    color: isActive ? '#4fc3f7' : 'rgba(255, 255, 255, 0.5)', 
                    // #4fc3f7 是亮藍色，未選中則是半透明白色 (變暗)

                    borderBottom: isActive ? '3px solid #4fc3f7' : '3px solid transparent', 
                    // 底部加一條線增強「選中感」

                    bgcolor: isActive ? 'rgba(79, 195, 247, 0.08)' : 'transparent',
                    // 選中時背景微微發光

                    '&:hover': {
                      color: '#fff', // 滑鼠移上去時變全白
                      bgcolor: 'rgba(255, 255, 255, 0.1)',
                    },
                  }}
                >
                  <Box sx={{ mb: 0.5 }}>{item.icon}</Box>
                  {item.label}
                </Button>
              );
            })}
          </Box>
        </Toolbar>
      </Container>
    </AppBar>
  );
};

export default Header;