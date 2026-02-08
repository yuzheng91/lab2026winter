import React from 'react';
import ReactDOM from 'react-dom/client';
import { BrowserRouter } from 'react-router-dom';
import { CssBaseline, ThemeProvider, createTheme } from '@mui/material'; // [新增 1] 引入 CssBaseline
import App from './App';
import './index.css';

// 建立一個預設主題 (可選，但建議加上)
const theme = createTheme({
  palette: {
    background: {
      default: "#f5f5f5" // 設定全站背景色，避免死白
    }
  }
});

ReactDOM.createRoot(document.getElementById('root')).render(
  <React.StrictMode>
    <BrowserRouter>
      <ThemeProvider theme={theme}>
        <CssBaseline /> {/* [新增 2] 這行就是消除白邊的關鍵！放在最外層 */}
        <App />
      </ThemeProvider>
    </BrowserRouter>
  </React.StrictMode>,
);