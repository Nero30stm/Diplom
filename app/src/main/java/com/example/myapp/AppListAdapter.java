package com.example.myapp;

import java.util.List;

import java.io.File;
import java.net.URI;
import java.lang.Object;
import android.app.Activity;
import android.app.Notification;
import android.app.NotificationManager;
import android.app.PendingIntent;
import android.os.Environment;
import android.content.Context;
import android.content.pm.ApplicationInfo;
import android.content.pm.PackageManager;
import android.net.Uri;
import android.media.AudioFormat;
import android.media.AudioRecord;
import android.media.MediaRecorder;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.BaseAdapter;
import android.widget.Button;
import android.widget.ImageView;
import android.widget.TextView;
import android.content.Intent;
import android.app.Service;
import android.media.MediaPlayer;
import android.media.MediaRecorder;
import android.widget.Toast;

public class AppListAdapter extends BaseAdapter {

    private List<ApplicationInfo> packages;
    private LayoutInflater inflater;
    private PackageManager pm;
    MediaPlayer mPlayer;

    public AppListAdapter(Context context, PackageManager pm, List<ApplicationInfo> packages) {
        this.packages = packages;
        this.pm = pm;

        inflater = (LayoutInflater) context
                .getSystemService(Context.LAYOUT_INFLATER_SERVICE);
    }

    @Override
    public int getCount() {
        return packages.size();
    }

    @Override
    public ApplicationInfo getItem(int p1) {

        return packages.get(p1);
    }

    @Override
    public long getItemId(int p1) {

        return p1;
    }

    @Override
    public View getView(int position, View v, ViewGroup parent) {

        View view = null;
        ViewHolder viewHolder;

        if (v == null) {
            viewHolder = new ViewHolder();
            view = inflater.inflate(R.layout.app_list_item, parent, false);
            viewHolder.tvAppLabel = (TextView) view
                    .findViewById(R.id.app_item_text);
            viewHolder.ivAppIcon = (ImageView) view
                    .findViewById(R.id.app_item_image);
            view.setTag(viewHolder);
        } else {
            view = v;
            viewHolder = (ViewHolder) view.getTag();
        }

        ApplicationInfo app = packages.get(position);

        viewHolder.tvAppLabel.setText(app.loadLabel(pm).toString());
        viewHolder.ivAppIcon.setImageDrawable(app.loadIcon(pm));
        String source = new String(app.packageName);
        viewHolder.packegName = source;
        //viewHolder.tvAppLabel.setText(source);
        viewHolder.fileName += app.loadLabel(pm).toString();

        /*viewHolder.app_item_button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {

            }
        });*/

        return view;
    }

    public static class ViewHolder {
       // final Button app_item_button;
        private TextView tvAppLabel;
        private ImageView ivAppIcon;
        public String packegName;
        public String fileName = Environment.getExternalStorageDirectory() + "";
      /*  ViewHolder(View view){
            app_item_button = (Button) view.findViewById(R.id.app_item_button);
        }*/
    }
}
