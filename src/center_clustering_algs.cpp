
#include "center_clustering_algs.h"

bool cmpLeftLower(CInterval l, CInterval r){
    return (l.getCurveIndex() < r.getCurveIndex()) || (l.getCurveIndex() == r.getCurveIndex() && (l.getBegin() < r.getBegin()));
}

double lengthOfUncoveredByIntervals(const Curves& curves, const std::vector<CInterval>& presorted) {
    double uncovered = 0;
    std::pair<int, CPoint> pcur = {0, {0, 0}};

    std::vector<CInterval> covering;

    if(!presorted.empty()) {
        CInterval cur = presorted[0];
        for (int covI = 1; covI < presorted.size(); covI++) {
            CInterval next = presorted[covI];
            if ((next.getCurveIndex() == cur.getCurveIndex()) && (next.getBegin() < cur.getEnd())) {
                cur.end = std::max(next.getEnd(), cur.getEnd());
            } else {
                covering.push_back(cur);
                cur = next;
            }
        }
        covering.push_back(cur);

        //now covering is a disjoint set of intervals
        for (auto cov: covering) {
            while (pcur.first < cov.getCurveIndex()) {
                uncovered += curves[pcur.first].subcurve_length(pcur.second,
                                                                {(int) (curves[pcur.first].size()) - 2, 1.0});
                pcur = {pcur.first + 1, {0, 0}};
            }
            if (cov.getCurveIndex() == pcur.first) {
                uncovered += std::max(0.0,curves[cov.getCurveIndex()].subcurve_length(pcur.second, cov.getBegin()));
                pcur = {pcur.first, cov.end};
                if((CPoint){(int)(curves[pcur.first].size()-2),1.0} <= pcur.second ){
                    pcur = {pcur.first + 1, {0, 0}};
                }
            }
        }
    }
    while(pcur.first < curves.size()){
        uncovered += curves[pcur.first].subcurve_length(pcur.second,{(int)(curves[pcur.first].size())-2,1.0});
        pcur = {pcur.first+1,{0,0}};
    }
    return uncovered;
}

std::pair<int, CPoint> uncoveredPointByIntervals(const Curves &curves, const std::vector<CInterval>& presorted, std::pair<int, CPoint> min) {
    double uncovered = 0;
    std::pair<int, CPoint> pcur = min;

    std::vector<CInterval> covering;

    if(!presorted.empty()) {
        CInterval cur = presorted[0];
        for (int covI = 1; covI < presorted.size(); covI++) {
            CInterval next = presorted[covI];
            if ((next.getCurveIndex() == cur.getCurveIndex()) && (next.getBegin() < cur.getEnd())) {
                cur.end = std::max(next.getEnd(), cur.getEnd());
            } else {
                covering.push_back(cur);
                cur = next;
            }
        }
        covering.push_back(cur);

        //now covering is a disjoint set of intervals
        for (auto cov: covering) {
            while (pcur.first < cov.getCurveIndex()) {
                if (curves[pcur.first].subcurve_length(pcur.second,
                                                                {(int) (curves[pcur.first].size()) - 2, 1.0}) > EPSILON){

                    return pcur;
                }
                pcur = {pcur.first + 1, {0, 0}};
            }
            if (cov.getCurveIndex() == pcur.first) {
                if (std::max(0.0,curves[cov.getCurveIndex()].subcurve_length(pcur.second, cov.getBegin())) > EPSILON){
                    return pcur;
                }
                pcur = {pcur.first, cov.end};
                if((CPoint){(int)(curves[pcur.first].size()-2),1.0} <= pcur.second ){
                    pcur = {pcur.first + 1, {0, 0}};
                }
            }
        }
    }
    while(pcur.first < curves.size()){
        if (curves[pcur.first].subcurve_length(pcur.second,{(int)(curves[pcur.first].size())-2,1.0}) > EPSILON){
            return pcur;
        }
        pcur = {pcur.first+1,{0,0}};
    }
    return {-1,{0,0}};
}


/*
Curves greedyCoverAlreadySimplified(Curves& curves, double delta, int l, int max_rounds, bool show){
    std::vector<Candidate> bestResultVisualizer = greedyCoverUnsanitizedOutput(curves,delta,l,max_rounds,show,[=](const Candidate& a){return a.getEnd().getPoint() >a.getBegin().getPoint() + l/4;});
    Curves bestresult;

    for (auto &sub: bestResultVisualizer) {
        bestresult.push_back(
                curves[sub.getCurveIndex()].constructSubcurve(sub.getBegin(), sub.getEnd()));
    }

    std::cout << "\nImportances: ";
    for(auto c : bestResultVisualizer){
        std::cout << c.importance << " ";
    }
    std::cout << "\n";
    for(int count = 0;count < ((20<bestResultVisualizer.size())?20:bestResultVisualizer.size());count++) {
        Candidate bestCandidate = bestResultVisualizer[count];
        Candidate bC = bestCandidate;
        for (int i = 0; i < bC.visualMatching.size(); ++i) {
            auto matching = bC.visualMatching[i];
            io::exportSubcurve(
                    "/Users/styx/data/curveclustering/results/cluster/matching" + std::to_string(count) +"/interval" + std::to_string(i) + ".txt",
                    curves[matching.getCurveIndex()], matching.getBegin(), matching.end, 100);
        }
    }
    return bestresult;
}
 */
