import React from 'react';
import { Routes, Route } from 'react-router-dom';
import { Box } from '@mui/material'; // 用來做整體佈局

// 引入元件
import Header from './components/Header';
import YeastAnalysis from './pages/YeastAnalysis';
import Wormbase from './pages/WormBase';
import { Home, Help, Contact } from './pages/StaticPages'; // 引入剛剛建立的靜態頁面

function App() {
  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', minHeight: '100vh' }}>
      {/* Header 放在 Routes 外面 
        這樣切換頁面時，Header 不會重新渲染，會一直停留在上面
      */}
      <Header />

      {/* 主內容區塊 
        flexGrow: 1 確保內容少時 footer (如果有) 會被推到底部，
        或是背景色能填滿剩餘空間
      */}
      <Box component="main" sx={{ flexGrow: 1, bgcolor: '#f5f5f5', minHeight: '100%' }}>
        <Routes>
          {/* 首頁 */}
          <Route path="/" element={<Home />} />
          
          {/* 功能頁 */}
          <Route path="/yeast" element={<YeastAnalysis />} />
          <Route path="/wormbase" element={<Wormbase />} />
          
          {/* 靜態頁 */}
          <Route path="/help" element={<Help />} />
          <Route path="/contact" element={<Contact />} />
          
          {/* 404 頁面 (捕捉所有未定義路由) */}
          <Route path="*" element={<Box sx={{p:4}}>404 Not Found</Box>} />
        </Routes>
      </Box>
    </Box>
  );
}

export default App;