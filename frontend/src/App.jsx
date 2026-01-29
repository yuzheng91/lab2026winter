import { Routes, Route } from 'react-router-dom';
import YeastAnalysis from './pages/YeastAnalysis';
import Wormbase from './pages/WormBase';

function App() {
  return (
    <Routes>
      {/* 當網址是 /yeast 時顯示分析頁面 */}
      <Route path="/yeast" element={<YeastAnalysis />} />
      <Route path="/wormbase" element={<Wormbase />} />
      {/* 你可以加入首頁或其他路由 */}
      <Route path="/" element={<div>Welcome to Gene Analysis System</div>} />
    </Routes>
  );
}

export default App;