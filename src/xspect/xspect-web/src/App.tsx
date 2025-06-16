import Header from "@/components/header"
import Landing from "@/components/landing";
import Classify from "@/components/classify"
import Filter from "@/components/filter"
import Models from "@/components/models"
import ResultPage from "@/components/result"
import FilteringResultPage from "@/components/filtering-result"
import { BrowserRouter as Router, Routes, Route } from "react-router-dom"
import ModelDetails from "./components/models-details";


function App() {
  return (
    <Router basename="/xspect-web">
      <div className="min-h-screen flex flex-col">
      <Header />
        <Routes>
          <Route path="/classify" element={<Classify />} />
          <Route path="/filter" element={<Filter />} />
          <Route path="/models" element={<Models />} />
          <Route path="/models/:model_slug" element={<ModelDetails />} />
          <Route path="/result/:classification_uuid" element={<ResultPage />} />
          <Route path="/filter-result/:filter_uuid" element={<FilteringResultPage />} />
          <Route path="*" element={<Landing />} />
        </Routes>
      </div>
    </Router>
  )
}

export default App
